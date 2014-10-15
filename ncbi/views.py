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


def get_csi_by_motif_type(curation_site_instances):
    """Return a tuple of motif-associated and non-motif-associated sites"""
    mas = curation_site_instances.filter(site_type="motif_associated")
    nmas = curation_site_instances.filter(site_type="non_motif_associated")
    return mas, nmas

def generate_tbl_string(curation_site_instances, test_export):
    """Given a collection of curation-site-instance objects, generate the .tbl
    text that will be sent to NCBI for RefSeq annotation. If test_export
    variable is True, the string is generated for test purposes only. No changes
    are done in the database."""
    tbl_str = ""
    tbl_str += write_header(curation_site_instances)
    # Create meta-sites
    mas, nmas = get_csi_by_motif_type(curation_site_instances)
    # Group sites by TF instances and create meta-site for each group,
    # individually.  This saves us to check TF instance of each
    # curation-site-instance object in the create_meta_site method.
    # Get all TF-instances/species values for all curation-site-instances.
    tf_instance_list = list(set(mas.values_list('curation__TF_instances')))
    print tf_instance_list
    # Find meta-sites for each group of curation-site-instances separately
    meta_sites = []
    for tf_insts in tf_instance_list:
        # Filter motif/non-motif associated sites specific to this TF instances
        filtered_m = mas.filter(curation__TF_instances=tf_insts)
        filtered_nm = nmas.filter(curation__TF_instances=tf_insts)
        meta_sites.extend(metasite.create_meta_sites(filtered_m, filtered_nm))
    assert meta_sites
    genome = meta_sites[0].genome
    # write each protein_bind feature
    for meta_site in meta_sites:
        # if any of the sites in meta-site is submitted to NCBI before, skip it
        ncbi_sites = [csi for csi in meta_site.cur_site_insts if
                      models.NCBISubmission.objects.filter(
                          genome_submitted_to=genome,
                          curation_site_instance=csi)]
        if ncbi_sites:
            continue
        create_ncbi_submission_object(meta_site, test_export)
        ncbi_sites.append(meta_site.delegate)
        tbl_str += write_location(meta_site)
        tbl_str += write_bound_moiety(meta_site)
        tbl_str += write_tf(meta_site)
        tbl_str += write_experimental_evidence(meta_site)
        tbl_str += write_regulation_evidence(meta_site)
        tbl_str += write_db_xref(meta_site)
    return str(tbl_str)

def write_header(curation_site_instances):
    """Return the header for the NCBI submission record."""
    return ('>Feature\tgi|%(gi)s|ref|%(accession)s\n' %
            {'gi': curation_site_instances[0].site_instance.genome.gi,
             'accession': curation_site_instances[0]\
             .site_instance.genome.genome_accession})

def create_ncbi_submission_object(meta_site, test_export):
    """If test_export is False, create a NCBI submission record in the database.
    If _test_export_ is True, do nothing, don't create anything in the database.
    """
    genome = meta_site.genome
    if not test_export:
        rec = models.NCBISubmission(genome_submitted_to=genome,
                                    curation_site_instance=meta_site.delegate)
        rec.save()

def write_location(meta_site):
    """Write binding site location"""
    start, end = (meta_site.delegate.site_instance.start+1,
                  meta_site.delegate.site_instance.end+1)
    if meta_site.delegate.site_instance.strand == -1:
        start, end = end, start
    return '%d\t%d\tprotein_bind' % (start, end) + '\n'

def write_bound_moiety(meta_site):
    """Write bound moiety to the NCBI record"""
    return '\t\t\tbound_moiety\t%s\n' % meta_site.delegate.curation.TF.name

def write_tf(meta_site):
    """Write TF name to the NCBI record."""
    # Check if submitted and reference organisms are same
    same_org = (meta_site.delegate.site_instance.genome.organism in
                [cur.site_species for cur in meta_site.curations])
    return ('\t\t\tnote\tTranscription factor binding site for %s. %s\n' %
            (meta_site.delegate.curation.TF_instances.all()[0].name,
             "" if same_org else
             "Originally reported in " + meta_site.curations[0].site_species))

def write_db_xref(meta_site):
    """Write the dbxref for the site"""
    return ('\t\t\tdb_xref\t%s\n' %
            dbxref.id2dbxref(int(meta_site.delegate.pk)))

def write_regulation_evidence(meta_site):
    """Return the string for NCBI submission containing evidence for gene
    regulation."""
    evidence_for_reg = set(reg.gene.locus_tag
                           for csi in meta_site.cur_site_insts
                           for reg in csi.regulation_set.all()
                           if reg.evidence_type == 'exp_verified')
    if evidence_for_reg:
        tbl_str = ('\t\t\tnote\tEvidence of regulation for: %s\n' %
                   (', '.join(evidence_for_reg)))
    else:
        tbl_str = ''
    return tbl_str

def write_experimental_evidence(meta_site):
    """Return the string for NCBI submission containing experimental evidence"""
    experiments = {}
    techniques = models.ExperimentalTechnique.objects\
                 .filter(preset_function__in=['binding', 'expression'])
    for exp in techniques:
        fcsis = [csi for csi in meta_site.cur_site_insts
                 if exp in csi.experimental_techniques.all()]
        experiments[exp] = list(set(csi.curation.publication.pmid
                                    for csi in fcsis))
    tbl_str = ''
    for exp, pmids in experiments.items():
        if not pmids: continue
        tbl_str += ('\t\t\texperiment\t%s [PMID: %s]\n' %
                    (exp.name, ', '.join(pmids)))
    return tbl_str

def download_full_tbl(genome_accession):
    """Download .tbl file for the given genome"""
    h = Entrez.efetch(db='nuccore', id=genome_accession,
                      rettype='ft', retmode='text')
    return h.read()

def download_full_fasta(genome_accession):
    """Download fasta file for the given genome"""
    h = Entrez.efetch(db='nuccore', id=genome_accession,
                      rettype='fasta', retmode='text')
    return h.read()

def download_genome_asn(genome_accession):
    """Download .asn file for the genome."""
    h = Entrez.efetch(db='nuccore', id=genome_accession,
                      rettype='null', retmode='text')
    return h.read()

def generate_readme_string():
    """Generate string to be written to readme file."""
    readme_str = ""
    readme_str += "Run: ./tbl2asn -p . -x .asn -R -Vvbr" + '\n'
    return readme_str

@user_passes_test(lambda u: u.is_staff)
def export_tbl_view(request):
    """Main view to export .tbl file for NCBI submission."""
    if request.method == "GET":
        form = ExportForm()
        return render(request, 'ncbi_export.html', {'form': form},
                      context_instance=RequestContext(request))

    # For given TF instance and genome, return all site instances in .tbl format
    form = ExportForm(request.POST)
    if not form.is_valid():
        return render(request, 'ncbi_export.html', {'form':form},
                      context_instance=RequestContext(request))

    genome_accession = form.cleaned_data['genome_accession']
    genome = models.Genome.objects.get(genome_accession=genome_accession)
    test_export = form.cleaned_data['is_test_export']
    return generate_zip_response(genome, test_export)

def generate_zip_response(genome, test_export):
    """Given a genome and the boolean whether if the export is for test
    purposes, generate a zip file containing asn file."""
    # get all curation_site_instances
    curation_site_instances = models.Curation_SiteInstance\
            .objects.filter(site_instance__genome=genome,
                            curation__NCBI_submission_ready=True,
                            site_type="motif_associated",
                            is_obsolete=False,
                            experimental_techniques__preset_function='binding')\
            .order_by('curation__TF_instances', 'site_instance__start')
    if len(curation_site_instances) == 0:
        msg = "No curation found for this genome."
        messages.add_message(request, messages.WARNING, msg)
        return render(request, 'ncbi_export.html', {'form':form},
                      context_instance=RequestContext(request))

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
