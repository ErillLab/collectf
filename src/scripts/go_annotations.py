"""Generate GO annotations.

TF-centric annotations
----------------------

For TF-binding annotations (enables) molecular function:

1. First case, one can simply specify that it is sequence specific:

- GO:0043565 sequence-specific DNA binding

which is what we would do if the paper simply reports binding and no evidence
of regulation of gene expression by the TF.

2. If there is evidence of TF-mediated regulation for at least one gene, but no
indication of whether the TF up- or down-regulates (repressor/activator) for
any gene in the curation, we can use the following term:

- GO:0001130 bacterial-type RNA polymerase transcription factor activity,
  sequence-specific DNA binding

- We also add explicitly: GO:0000976 transcription regulatory region
  sequence-specific DNA binding

meaning that we tag along a GO:0000976 to the other terms (i.e. 2 annotations
generated)

3. And if there is evidence of regulation AND at least one gene is shown to be
activated (TF: activator) or repressed (TF: repressor), we can use:

- GO:0001216 bacterial-type RNA polymerase transcriptional activator activity,
  sequence-specific DNA binding

- GO:0001217 bacterial-type RNA polymerase transcriptional repressor activity,
  sequence-specific DNA binding

- we also add explicitly GO:0000976 transcription regulatory region
  sequence-specific DNA binding to both of the above

meaning that we tag along a GO:0000976 to the other terms (i.e. 2 annotations
generated)

4. For regulation annotations (involved_in) biological process:

Here, for a reported but unspecified regulatory effect (e.g. no indication of
activation vs. repression), we have:

- GO:0006355 regulation of transcription, DNA-templated

and for activation/repression, we get:

- GO:0045893 positive regulation of transcription, DNA-templated

- GO:0045893 positive regulation of transcription, DNA-templated

Again, this means two annotations if the TF activates at least one gene and
represses at least one gene.
"""

import os
import time

from django.conf import settings

from core import models
from core import entrez_utils


def new_TF_centric_gpad_entry(uniprot_id, go_id, db_ref,
                              evidence_code, relationship):
    """Returns a dictionary for a TF-centric GO annotation.

    See http://geneontology.org/page/gene-product-association-data-gpad-format
    """
    return dict(db='UniProtKB',
                db_object_id=uniprot_id,
                relationship=relationship,
                go_id=go_id,
                db_ref=('PMID:' + db_ref),
                evidence_code=evidence_code,
                with_from='',
                interacting_taxon_id='',
                date=time.strftime('%Y%m%d'),
                assigned_by='CollecTF',
                annotation_extension='',
                annotation_properties='')


def new_gene_centric_gpad_entry(refseq_acc, go_id, db_ref, evidence_code):
    """Returns a dictionary for a gene-centric GO annotation.

    See http://geneontology.org/page/gene-product-association-data-gpad-format
    """
    return dict(db='RefSeq',
                db_object_id=refseq_acc,
                relationship='involved_in',
                go_id=go_id,
                db_ref=('PMID:' + db_ref),
                evidence_code=evidence_code,
                with_from='',
                interacting_taxon_id='',
                date=time.strftime('%Y%m%d'),
                assigned_by='CollecTF',
                annotation_extension='',
                annotation_properties='')


def get_exp_regulations(TF_instance):
    """Returns the experimentally-validated regulations of a given TF."""
    return models.Regulation.objects.filter(
        evidence_type='exp_verified',
        curation_site_instance__curation__TF_instances=TF_instance)


def get_inferred_regulations(TF_instance):
    """Returns inferred regulations of a given TF."""
    return models.Regulation.objects.filter(
        evidence_type='inferred',
        curation_site_instance__curation__TF_instances=TF_instance)


def TF_centric_binding_annotations(TF_instance):
    """Generates list of binding GO annotations for the given TF instance."""
    print TF_instance
    uniprot_id = TF_instance.uniprot_accession
    go_annotations = []
    annotated_pmids = set()
    exp_regulations = get_exp_regulations(TF_instance)
    inferred_regulations = get_inferred_regulations(TF_instance)

    # If there is evidence of regulation AND at least one gene is shown to
    # be activated or repressed,
    for reg in (r for r in exp_regulations if r.mode in ['ACT', 'REP']):
        pmid = reg.ref_pmid
        for tech in reg.binding_experimental_techniques:
            if reg.mode == 'ACT':
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0001216', pmid, tech.EO_term, 'enables'))
            else:  # reg.mode = 'REP'
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0001217', pmid, tech.EO_term, 'enables'))
            go_annotations.append(new_TF_centric_gpad_entry(
                uniprot_id, 'GO:0000976', pmid, tech.EO_term, 'enables'))
        annotated_pmids.add(pmid)

    # If there is evidence of TF-mediated regulation for at least one gene,
    # but no indication of whether the TF up- or down-regulates
    # (repressor/activator) for any gene in the curation
    for reg in (r for r in exp_regulations if r.mode not in ['ACT', 'REP']):
        pmid = reg.ref_pmid
        if pmid not in annotated_pmids:
            for tech in reg.binding_experimental_techniques:
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0001130', pmid, tech.EO_term, 'enables'))
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0000976', pmid, tech.EO_term, 'enables'))
                annotated_pmids.add(pmid)

    # If there is no regulation evidence, simply report binding
    for reg in inferred_regulations:
        pmid = reg.ref_pmid
        if pmid not in annotated_pmids:
            for tech in reg.binding_experimental_techniques:
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0043565', pmid, tech.EO_term, 'enables'))
                annotated_pmids.add(pmid)

    return set(map(gpad_entry_to_str, go_annotations))


def TF_centric_regulation_annotations(TF_instance):
    """Generates list of regulation GO annotations for the given TF instance"""
    uniprot_id = TF_instance.uniprot_accession
    go_annotations = []
    for reg in get_exp_regulations(TF_instance):
        for tech in reg.expression_experimental_techniques:
            if reg.mode == 'ACT':
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0045893', reg.ref_pmid, tech.EO_term,
                    'involved_in'))
            elif reg.mode == 'REP':
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0045892', reg.ref_pmid, tech.EO_term,
                    'involved_in'))
            else:
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0006355', reg.ref_pmid, tech.EO_term,
                    'involved_in'))
    return set(map(gpad_entry_to_str, go_annotations))


def batch_refseq_accession(genome_accession):
    """Retrieves product NP/YP/WP accession numbers of the given genome.

    Returns dictionary {locus_tag: RefSeq_acc}.
    """
    print "Downloading", genome_accession
    genome_rec = entrez_utils.get_genome(genome_accession)
    genes = entrez_utils.get_genes(genome_rec)
    gene_dict = {}
    for gene in genes:
        if gene['protein_id'] and len(gene['protein_id']) == 1:
            protein_id = gene['protein_id'][0].split('.')[0]
            gene_dict[gene['locus_tag']] = protein_id
            for old_locus_tag in gene['old_locus_tag']:
                gene_dict[old_locus_tag] = protein_id
    return gene_dict


def gene_centric_annotations(TF_instance):
    if not TF_instance.GO_term:
        return []
    print TF_instance
    go_annotations = []
    exp_regulations = get_exp_regulations(TF_instance)
    regulated_genomes = set(reg.gene.genome for reg in exp_regulations
                            if reg.expression_experimental_techniques)
    for genome in regulated_genomes:
        regulations_in_genome = exp_regulations.filter(gene__genome=genome)
        refseq_accs = batch_refseq_accession(genome.genome_accession)
        for reg in regulations_in_genome:
            locus_tag = reg.gene.locus_tag
            print locus_tag
            for tech in reg.expression_experimental_techniques:
                if refseq_accs.get(locus_tag):
                    go_annotations.append(new_gene_centric_gpad_entry(
                        refseq_accs[locus_tag], TF_instance.GO_term.GO_term_id,
                        reg.ref_pmid, tech.EO_term))
    print len(go_annotations)
    return set(map(gpad_entry_to_str, go_annotations))


def gpad_entry_to_str(gpad_entry):
    return '\t'.join(gpad_entry[field] for field in
                     ['db', 'db_object_id', 'relationship', 'go_id', 'db_ref',
                      'evidence_code', 'with_from', 'interacting_taxon_id',
                      'date', 'assigned_by', 'annotation_extension',
                      'annotation_properties'])


def generate_gpad_file():
    """Collects GO annotations for all TF-instances and generates GPAD file."""
    header = [
        "!gpa-version: 1.1",
        "!Submission Date: " + time.strftime('%m/%d/%Y'),
        "!",
        "!Project_name: CollecTF - Experimentally validated",
        "!URL: http://www.collectf.org/",
        "!Contact Email: collectf@umbc.edu",
        "!",
        '\t'.join(['!DB', 'DB_Object_ID', 'Relationship', 'GO ID', 'Reference',
                   'Evidence code', 'With (or) From', 'Interacting taxon ID',
                   'Date', 'Assigned By', 'Annotation Extension',
                   'Annotation Properties'])]
    content = []
    for TF_instance in models.TFInstance.objects.all():
        content.extend(TF_centric_binding_annotations(TF_instance))
        content.extend(TF_centric_regulation_annotations(TF_instance))
        content.extend(gene_centric_annotations(TF_instance))
    export_file = os.path.join(settings.STATICFILES_DIRS[0], 'collectf.gpad')
    with open(export_file, 'w') as f:
        f.write('\n'.join(header + content))


def run():
    generate_gpad_file()
