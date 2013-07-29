# Create your views here.
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render
from django.shortcuts import render_to_response
from django.template import RequestContext
import models
import os
from forms import ExportForm
from baseapp import utils
from baseapp import bioutils
import tempfile

@user_passes_test(lambda u: u.is_staff)
def export_tbl_view(request):
    if request.method=="GET":
        form = ExportForm()
        return render(request,
                      'ncbi_export.html',
                      {'form': form},
                      context_instance=RequestContext(request))

    else:
        form = ExportForm(request.POST)
        form.is_valid()
        tbl_data = export_tbl(form.cleaned_data['TF_instances'], form.cleaned_data['genomes'])
        return render_to_response('ncbi_export_result.html',
                                  {
                                      'tbl': tbl_data,
                                      'TF_instance': form.cleaned_data['TF_instances'],
                                      'genome': form.cleaned_data['genomes'],
                                  },
                                  context_instance=RequestContext(request))


def export_tbl(TF_instance, genome):
    """For given TF instance and genome, return all site instances in .tbl format"""
    # open tbl file to write
    tbl_file = tempfile.TemporaryFile()
    # get all curation_site_instances
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome=genome,
        curation__TF_instance=TF_instance,
        curation__NCBI_submission_ready=True,
        is_motif_associated=True,
        curation__experimental_techniques__preset_function__in=['binding', 'expression'])
    
    if len(set(csi.curation.TF for csi in curation_site_instances)) > 1:
        msg = "Inconsistent TF-TF_instance links. This TF_instance is linked to more than one TF."
        tbl_file.write(msg + '\n')
        tbl_file.seek(0) # goto beginnning of the file
        return tbl_file.read()

    meta_sites = dict()
    for csi in curation_site_instances:
        # search for a meta-site-instance
        for i in meta_sites.keys():
            if bioutils.overlap_site_meta_site(csi.site_instance, [m.site_instance for m in meta_sites[i]]):
                meta_sites[i].append(csi)
                break
        else:
            meta_sites[len(meta_sites)+1] = [csi]

    # FIXME
    tbl_file.write('>Feature prot_%s_genome_%s\n' % (TF_instance.protein_accession, genome.genome_accession))
    for ms_id, meta_site in meta_sites.items():
        # pick a site to report its id as dbxref.
        ncbi_sites = [csi for csi in meta_site if csi.is_ncbi_xref]
        if not ncbi_sites: # pick first site as ncbi_xref
            ncbi_sites.append(meta_site[0])
        ncbi_site = ncbi_sites[0]
        tbl_file.write('%d\t%d\tprotein_bind' % (ncbi_site.site_instance.start, ncbi_site.site_instance.end) + '\n')
        tbl_file.write('\t\t\tbound_moiety\t%s\n' % ncbi_site.curation.TF.name)
        tbl_file.write('\t\t\tnote\tTranscription factor binding site\n')
        # write experimental evidences
        experiments = {}
        for exp in models.ExperimentalTechnique.objects.filter(preset_function__in=['binding', 'expression']):
            filtered_csis = [csi for csi in meta_site if exp in csi.curation.experimental_techniques.all()]
            experiments[exp] = list(set([csi.curation.publication.pmid for csi in filtered_csis]))
        for exp,pmids in experiments.items():
            if not pmids: continue
            tbl_file.write('\t\t\texperiment\t%s [PMID: %s]\n' % (exp.name, ', '.join(pmids)))

        # write regulation note
        evidence4regulation = set([reg.gene.locus_tag for csi in meta_site
                                   for reg in csi.regulation_set.all()
                                   if reg.evidence_type=="exp_verified"])
        if evidence4regulation:
            tbl_file.write('\t\t\tnote\tEvidence of regulation for: %s\n' % (', '.join(evidence4regulation)))

        # write dbxref
        tbl_file.write('\t\t\tdb_xref\t%s\n' % utils.id2dbxref(int(meta_site[0].pk)))

    tbl_file.seek(0) # goto beginnning of the file
    return tbl_file.read()


