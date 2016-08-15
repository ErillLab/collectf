"""Script to generate NCBI genome annotations."""

import os


from core import models
from core import metasite
from core import dbxref

from django.conf import settings

from tqdm import tqdm


def get_site_instances(accession_number):
    """Returns the meta-sites, given the genome accession number."""
    genome = models.Genome.objects.get(genome_accession=accession_number)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome=genome,
        curation__NCBI_submission_ready=True,
        site_type='motif_associated',
        is_obsolete=False,
        experimental_techniques__preset_function='binding').order_by(
            'site_instance__start').distinct()

    # Get list of TF instances
    TF_instances = list(set(curation_site_instances.values_list(
        'curation__TF_instances')))
    # Create metasites
    meta_sites = []
    for TF_instance in TF_instances:
        meta_sites.extend(metasite.create_meta_sites(
            curation_site_instances.filter(curation__TF_instances=TF_instance)))

    gff_lines = [create_gff_line(meta_site) for meta_site in meta_sites]
    return '\n'.join(gff_lines)


def create_gff_line(meta_site):
    """Generates a GFF line from the given meta-site."""
    def attr_bound_moiety():
        return 'bound_moiety: %s' % meta_site.delegate.curation.TF.name

    def attr_TF():
        return 'TF: %s' % meta_site.delegate.curation.TF_instances.all()[0].refseq_accession

    def attr_dbxref():
        return 'db_xref: %s' % dbxref.to_ncbi_dbxref(meta_site.delegate.pk)

    def attr_regulation_evidence():
        regulated_genes = set(
            reg.gene.locus_tag
            for csi in meta_site.curation_site_instances
            for reg in csi.regulation_set.all()
            if reg.evidence_type == 'exp_verified')
        return ('Evidence of regulation for %s.' % ', '.join(regulated_genes)
                if regulated_genes else '-')

    def attr_experimental_evidence():
        experiments = {}
        techniques = models.ExperimentalTechnique.objects.filter(
            preset_function__in=['binding', 'expression'])
        for exp in techniques:
            cur_site_insts = [csi for csi in meta_site.curation_site_instances
                              if exp in csi.experimental_techniques.all()]
            pmids = set(csi.curation.publication.pmid for
                        csi in cur_site_insts)
            if pmids:
                experiments[exp] = pmids

        evidence_str = '; '.join(('experiment: %s (%s) (PMID: %s)' %
                                  (exp.name, exp.EO_term, ', '.join(pmids)))
                                 for exp, pmids in experiments.items())
        return evidence_str

    delegate_site_instance = meta_site.delegate_site_instance

    fields = [delegate_site_instance.genome.genome_accession,
              'CollecTF',
              'protein_bind',
              str(delegate_site_instance.start),
              str(delegate_site_instance.end),
              '',               # Empty score
              '+' if delegate_site_instance.strand == 1 else '-',
              '',               # Frame not applicable
              '; '.join([attr_bound_moiety(),
                         attr_TF(),
                         attr_dbxref(),
                         attr_regulation_evidence(),
                         attr_experimental_evidence()])]
    return '\t'.join(fields)


def list_of_genomes():
    """Returns the list of genome accession numbers that CollecTF has data for."""
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__NCBI_submission_ready=True,
        site_type='motif_associated',
        is_obsolete=False,
        experimental_techniques__preset_function='binding')
    genomes = curation_site_instances.values_list(
        'site_instance__genome__genome_accession', flat=True).distinct()
    # Export the list of genome accession numbers.
    export_file = os.path.join(settings.STATICFILES_DIRS[0], 'ncbi', 'accessions.txt')
    with open(export_file, 'w') as f:
        f.write('\n'.join(genomes))

    return genomes


def run():
    genome_list = list_of_genomes()
    for genome in tqdm(genome_list):
        gff = get_site_instances(genome)
        export_file = os.path.join(settings.STATICFILES_DIRS[0], 'ncbi', genome+'.gff')
        with open(export_file, 'w') as f:
            f.write(gff)
