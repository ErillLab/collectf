# Create your views here.
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponse
import models
import os
from StringIO import StringIO
from zipfile import ZipFile
from forms import ExportForm
from baseapp import utils
from baseapp import bioutils

def generate_tbl_string(curation_site_instances):
    tbl_str = ""
    # Create meta-sites
    meta_sites = dict()
    for csi in curation_site_instances:
        # search for a meta-site-instance
        for i in meta_sites.keys():
            if bioutils.overlap_site_meta_site(csi, meta_sites[i]):
                meta_sites[i].append(csi)
                break
        else:
            meta_sites[len(meta_sites)+1] = [csi]
    # FIXME
    # header
    tbl_str += ('>Feature\tgi|%(gi)s|ref|%(accession)s|\tTF_%(TF)s_%(TF_accession)s\n' %
                {'gi': curation_site_instances[0].site_instance.genome.gi,
                 'accession': curation_site_instances[0].site_instance.genome.genome_accession,
                 'TF': curation_site_instances[0].curation.TF.name,
                 'TF_accession': curation_site_instances[0].curation.TF_instance.protein_accession})
    # write each protein_bind feature
    for ms_id, meta_site in meta_sites.items():
        # pick a site to report its id as dbxref.
        ncbi_sites = [csi for csi in meta_site if csi.ncbi_submission]
        if not ncbi_sites: # pick first site as ncbi_xref
            n = models.NCBISubmission(is_obsolete=False, why_obsolete="")
            n.save()
            meta_site[0].ncbi_submission = n
            meta_site[0].save()
            ncbi_sites.append(meta_site[0])
        ncbi_site = ncbi_sites[0]
        tbl_str += ('%d\t%d\tprotein_bind' % (ncbi_site.site_instance.start, ncbi_site.site_instance.end) + '\n')
        tbl_str += ('\t\t\tbound_moiety\t%s\n' % ncbi_site.curation.TF.name)
        tbl_str += ('\t\t\tnote\tTranscription factor binding site\n')
        # write experimental evidences
        experiments = {}
        for exp in models.ExperimentalTechnique.objects.filter(preset_function__in=['binding', 'expression']):
            filtered_csis = [csi for csi in meta_site if exp in csi.curation.experimental_techniques.all()]
            experiments[exp] = list(set([csi.curation.publication.pmid for csi in filtered_csis]))
        for exp,pmids in experiments.items():
            if not pmids: continue
        # write regulation note
        evidence4regulation = set([reg.gene.locus_tag for csi in meta_site
                                   for reg in csi.regulation_set.all()
                                   if reg.evidence_type=="exp_verified"])
        if evidence4regulation:
            tbl_str += ('\t\t\tnote\tEvidence of regulation for: %s\n' % (', '.join(evidence4regulation)))
        # write dbxref
        tbl_str += ('\t\t\tdb_xref\t%s\n' % utils.id2dbxref(int(meta_site[0].pk)))

    return str(tbl_str)

def generate_src_string(curation_site_instances):
    src_str = ""
    src_str += ("sequence_id\torganism\tchromosome\n")
    src_str += ("gi|%(gi)s|ref|%(accession)s|\t%(organism)s\t%(chromosome)s\n" %
                {'gi': curation_site_instances[0].site_instance.genome.gi,
                 'accession': curation_site_instances[0].site_instance.genome.genome_accession,
                 'organism': curation_site_instances[0].site_instance.genome.organism,
                 'chromosome': curation_site_instances[0].site_instance.genome.chromosome})
    return str(src_str)

def generate_readme_string():
    readme_str = ""
    readme_str += "Create a submission template at http://www.ncbi.nlm.nih.gov/WebSub/template.cgi\n"
    readme_str += "Gneome file (.fsa) and feature table (.tbl) must have the same filename prefix\n"
    readme_str += "Run: ./tbl2asn -p. -t <template_file> -i <genome FASTA file> -V vbr\n"
    return readme_str

@user_passes_test(lambda u: u.is_staff)
def export_tbl_view(request):
    if request.method=="GET":
        form = ExportForm()
        return render(request,
                      'ncbi_export.html',
                      {'form': form},
                      context_instance=RequestContext(request))

    # For given TF instance and genome, return all site instances in .tbl format
    form = ExportForm(request.POST)
    form.is_valid()
    TF_instance = form.cleaned_data['TF_instances']
    genome = form.cleaned_data['genomes']

    # get all curation_site_instances
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome=genome,
        curation__TF_instance=TF_instance,
        curation__NCBI_submission_ready=True,
        is_motif_associated=True,
        curation__experimental_techniques__preset_function__in=['binding', 'expression'])

    if len(set(csi.curation.TF for csi in curation_site_instances)) > 1:
        msg = "Inconsistent TF-TF_instance links. This TF_instance is linked to more than one TF."
        return render(request, 'ncbi_export.html', {'form':form}, context_instance=RequestContext(request))

    if len(curation_site_instances) == 0:
        msg = "No curation found for this TF and genome."
        return render(request, 'ncbi_export.html', {'form':form}, context_instance=RequestContext(request))

    tbl_str = generate_tbl_string(curation_site_instances)
    src_str = generate_src_string(curation_site_instances)
    readme_str = generate_readme_string()

    # create a zip file
    filename = 'testing'
    in_memory = StringIO()
    zip = ZipFile(in_memory, 'a')
    zip.writestr(filename+'.tbl', tbl_str)
    zip.writestr(filename+'.src', src_str)
    zip.writestr("README.txt", readme_str)
    zip.close()
    response = HttpResponse(content_type='application/zip')
    response['Content-Disposition'] = 'attachment;filename=tbl_export.zip'
    in_memory.seek(0)
    response.write(in_memory.read())
    return response


