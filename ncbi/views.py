from django.shortcuts import render
from forms import ExportForm
from base import models
from django.contrib.auth.decorators import user_passes_test
from django.template import RequestContext
from django.contrib import messages
from django.http import HttpResponse
from base import metasite
from base.templatetags import dbxref
from StringIO import StringIO
from zipfile import ZipFile
from Bio import Entrez
Entrez.email = 'sefa1@umbc.edu'
import template_sbt


def generate_tbl_string(curation_site_instances, test_export):
    """Given a collection of curation-site-instance objects, generate the .tbl
    text that will be sent to NCBI for RefSeq annotation. If test_export
    variable is True, the string is generated for test purposes only. No changes
    are done in the database."""
    tbl_str = ""
    # Create meta-sites
    motif_associated = curation_site_instances.filter(site_type="motif_associated")
    non_motif_associated = curation_site_instances.filter(site_type="non_motif_associated")
    meta_sites = metasite.create_meta_sites(motif_associated, non_motif_associated)

    # header
    tbl_str += ('>Feature\tgi|%(gi)s|ref|%(accession)s\n' %
                {'gi': curation_site_instances[0].site_instance.genome.gi,
                 'accession': curation_site_instances[0].site_instance.genome.genome_accession})

    assert meta_sites
    genome = meta_sites[0].genome
    # write each protein_bind feature
    for meta_site in meta_sites:
        # if any of the sites in the meta-site is submitted to NCBI before, skip it.
        ncbi_sites = [csi for csi in meta_site.cur_site_insts
                      if models.NCBISubmission.objects.filter(
                              genome_submitted_to=genome,
                              curation_site_instance=csi)]
        if ncbi_sites:
            continue

        # Pick up the first site as ncbi_Xref
        if not test_export:
            n = models.NCBISubmission(genome_submitted_to=genome, curation_site_instance=meta_site.delegate)
            n.save()

        ncbi_sites.append(meta_site.delegate)
        start, end = meta_site.delegate.site_instance.start+1, meta_site.delegate.site_instance.end+1
        if meta_site.delegate.site_instance.strand == -1:
            start,end = end,start

        tbl_str += ('%d\t%d\tprotein_bind' % (start,end) + '\n')
        tbl_str += ('\t\t\tbound_moiety\t%s\n' % meta_site.delegate.curation.TF.name)
        tbl_str += ('\t\t\tnote\tTranscription factor binding site for %s\n' % meta_site.delegate.curation.TF_instances.all()[0].name)
        # write experimental evidences
        experiments = {}
        for exp in models.ExperimentalTechnique.objects.filter(preset_function__in=['binding', 'expression']):
            filtered_csis = [csi for csi in meta_site.cur_site_insts if exp in csi.experimental_techniques.all()]
            experiments[exp] = list(set([csi.curation.publication.pmid for csi in filtered_csis]))
        for exp,pmids in experiments.items():
            if not pmids: continue
            tbl_str += ('\t\t\texperiment\t%s [PMID: %s]\n' % (exp.name, ', '.join(pmids)))
        # write regulation note
        evidence4regulation = set([reg.gene.locus_tag for csi in meta_site.cur_site_insts
                                   for reg in csi.regulation_set.all()
                                   if reg.evidence_type=="exp_verified"])
        if evidence4regulation:
            tbl_str += ('\t\t\tnote\tEvidence of regulation for: %s\n' % (', '.join(evidence4regulation)))
        # write dbxref
        tbl_str += ('\t\t\tdb_xref\t%s\n' % dbxref.id2dbxref(int(meta_site.delegate.pk)))
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

def download_full_tbl(genome_accession):
    h = Entrez.efetch(db='nuccore', id=genome_accession, rettype='ft', retmode='text')
    return h.read()

def download_full_fasta(genome_accession):
    h = Entrez.efetch(db='nuccore', id=genome_accession, rettype='fasta', retmode='text')
    return h.read()

def download_genome_asn(genome_accession):
    h = Entrez.efetch(db='nuccore', id=genome_accession, rettype='null', retmode='text')
    return h.read()

def generate_readme_string():
    readme_str = ""
    readme_str += "Run: ./tbl2asn -p . -x .asn -R -Vvbr" + '\n'
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
    genome = form.cleaned_data['genomes']
    test_export = form.cleaned_data['is_test_export']

    # get all curation_site_instances
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome=genome,
        curation__NCBI_submission_ready=True,
        site_type="motif_associated",
        is_obsolete=False,
        experimental_techniques__preset_function='binding').order_by('curation__TF_instances',
                                                                     'site_instance__start')

    if len(curation_site_instances) == 0:
        msg = "No curation found for this genome."
        messages.add_message(request, messages.WARNING, msg)
        return render(request, 'ncbi_export.html', {'form':form}, context_instance=RequestContext(request))

    collectf_tbl_str = generate_tbl_string(curation_site_instances, test_export)
    genome_asn_str = download_genome_asn(genome.genome_accession)
    readme_str = generate_readme_string()

    # create a zip file
    filename = genome.genome_accession.replace('.', '_')
    in_memory = StringIO()
    zip = ZipFile(in_memory, 'a')
    zip.writestr(filename+'.tbl', collectf_tbl_str)
    zip.writestr(filename+'.asn', genome_asn_str)
    zip.writestr("CollecTF_template.sbt", template_sbt.COLLECTF_TEMPLATE_SBT)
    zip.writestr("command.txt", readme_str)
    zip.close()
    response = HttpResponse(content_type='application/zip')
    response['Content-Disposition'] = 'attachment;filename=%s.zip' % filename
    in_memory.seek(0)
    response.write(in_memory.read())
    return response
