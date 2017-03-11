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

03/09/2017
- Modified script so that it adds WITH field for IPI-derived codes, using the
  genome accession from RefSeq as ID for WITH field
- Modified script so that it also generates a cellular component annotation
  using GO:0032993 - protein-DNA complex for binding-supported experiments.
"""

import os
import time

from django.conf import settings

import uniprot

from core import models
from core import entrez_utils


def new_TF_centric_gpad_entry(uniprot_id, go_id, db_ref,
                              evidence_code, relationship, with_field):
    """Returns a dictionary for a TF-centric GO annotation.

    See http://geneontology.org/page/gene-product-association-data-gpad-format
    """
    return dict(db='UniProtKB',
                db_object_id=uniprot_id,
                relationship=relationship,
                go_id=go_id,
                db_ref=('PMID:' + db_ref),
                evidence_code=evidence_code,
                with_from=with_field,
                interacting_taxon_id='',
                date=time.strftime('%Y%m%d'),
                assigned_by='CollecTF',
                annotation_extension='',
                annotation_properties='')


def new_gene_centric_gpad_entry(uniprot_id, go_id, db_ref, evidence_code, with_field):
    """Returns a dictionary for a gene-centric GO annotation.

    See http://geneontology.org/page/gene-product-association-data-gpad-format
    """
    return dict(db='UniProtKB',
                db_object_id=uniprot_id,
                relationship='involved_in',
                go_id=go_id,
                db_ref=('PMID:' + db_ref),
                evidence_code=evidence_code,
                with_from=with_field,
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

def get_ECO_IPI_terms():
    """Returns list of ECO terms belonging to IPI, which require
       WITH field (with genome ID) in the GO annotation"""
    return ['ECO:0005618', 'ECO:0005621', 'ECO:0005626', 'ECO:0005630', \
            'ECO:0005631', 'ECO:0005635', 'ECO:0005643', 'ECO:0005647', \
            'ECO:0005656', 'ECO:0005665']

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
            with_field='' #define WITH field, fill it up only if IPI evidence used
            if tech.EO_term in get_ECO_IPI_terms():
                with_field='RefSeq:'+reg.curation_site_instance.site_instance.genome.genome_accession
            if reg.mode == 'ACT':
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0001216', pmid, tech.EO_term, 'enables',with_field))
            else:  # reg.mode = 'REP'
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0001217', pmid, tech.EO_term, 'enables',with_field))
            go_annotations.append(new_TF_centric_gpad_entry(
                uniprot_id, 'GO:0000976', pmid, tech.EO_term, 'enables',with_field))
            go_annotations.append(new_TF_centric_gpad_entry(
                uniprot_id, 'GO:0032993', pmid, tech.EO_term, 'part_of',with_field))
        annotated_pmids.add(pmid)

    # If there is evidence of TF-mediated regulation for at least one gene,
    # but no indication of whether the TF up- or down-regulates
    # (repressor/activator) for any gene in the curation
    for reg in (r for r in exp_regulations if r.mode not in ['ACT', 'REP']):
        pmid = reg.ref_pmid
        if pmid not in annotated_pmids:
            for tech in reg.binding_experimental_techniques:
                with_field='' #define WITH field, fill it up only if IPI evidence used
                if tech.EO_term in get_ECO_IPI_terms():
                    with_field='RefSeq:'+reg.curation_site_instance.site_instance.genome.genome_accession            
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0001130', pmid, tech.EO_term, 'enables',with_field))
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0000976', pmid, tech.EO_term, 'enables',with_field))
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0032993', pmid, tech.EO_term, 'part_of',with_field))
                annotated_pmids.add(pmid)

    # If there is no regulation evidence, simply report binding
    for reg in inferred_regulations:
        pmid = reg.ref_pmid
        if pmid not in annotated_pmids:
            for tech in reg.binding_experimental_techniques:
                with_field='' #define WITH field, fill it up only if IPI evidence used
                if tech.EO_term in get_ECO_IPI_terms():
                    with_field='RefSeq:'+reg.curation_site_instance.site_instance.genome.genome_accession                
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0043565', pmid, tech.EO_term, 'enables',with_field))
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0032993', pmid, tech.EO_term, 'part_of',with_field))
                annotated_pmids.add(pmid)

    return set(map(gpad_entry_to_str, go_annotations))


def TF_centric_regulation_annotations(TF_instance):
    """Generates list of regulation GO annotations for the given TF instance"""
    uniprot_id = TF_instance.uniprot_accession
    go_annotations = []
    for reg in get_exp_regulations(TF_instance):
        for tech in reg.expression_experimental_techniques:
            with_field='' #define WITH field, fill it up only if IPI evidence used
            if tech.EO_term in get_ECO_IPI_terms():
                with_field='RefSeq:'+reg.curation_site_instance.site_instance.genome.genome_accession            
            if reg.mode == 'ACT':
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0045893', reg.ref_pmid, tech.EO_term,
                    'involved_in',with_field))
            elif reg.mode == 'REP':
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0045892', reg.ref_pmid, tech.EO_term,
                    'involved_in',with_field))
            else:
                go_annotations.append(new_TF_centric_gpad_entry(
                    uniprot_id, 'GO:0006355', reg.ref_pmid, tech.EO_term,
                    'involved_in',with_field))
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
            protein_id = gene['protein_id'][0]#.split('.')[0]
            gene_dict[gene['locus_tag']] = protein_id
            for old_locus_tag in gene['old_locus_tag']:
                gene_dict[old_locus_tag] = protein_id
    return gene_dict

def uniprot_to_refseq(genome_acc):
    """Gets a genome accession. Queries uniprot to get all UniProt accessions
       associated with genome. Downloads genome from NCBI and gets all refseq 
       protein records associated with it.
       For each uniprot record, checks whether it maps to a refseq record.
       If so, it returns a dictionary with the complete mapping.
    """
    
    #make uniprot query to get all uniprot IDs for a given NCBI genome accession
    #the result is something like: [u'A0A0H3IFM0', u'B8H3X3', u'A0A0H3IXV8']
    u_accs=uniprot.map(genome_acc, f='REFSEQ_NT_ID' , t='ACC')
    
    if len(u_accs)>0:
        uniprot_accs = list(u_accs.get(genome_acc,[]))
    
    #download genome and get refseq protein accessions
    #the result is something like: {'CCNA_00250': 'YP_002515625', 'CCNA_02456': 'YP_002517829'}
    refseq_accs=batch_refseq_accession(genome_acc)
    
    gene_prot_dict = {}
    #get mapping from all uniprot records (associated to genome) to refseq ids
    #the result is something like:
    #{u'A0A0H3IFM0': {u'WP_024265740.1', u'YP_008877611.1'},
    # u'A0A0H3IXV8': {u'WP_024265959.1', u'YP_009020565.1'}}
    if len(uniprot_accs)>0:
        uni_refseqs=uniprot.map(uniprot_accs, f='ACC', t='P_REFSEQ_AC')
            
        #for each gene with mapped refseq protein record
        for locus_tag, refseq_acc in refseq_accs.iteritems():
            #check whether the refseq protein record has also been mapped from uniprot
            #for each uniprot record mapped to a refseq record
            for uni, ref in uni_refseqs.iteritems():
                if refseq_acc in ref:
                    protein_ids={'refseq':refseq_acc, 'uniprot':uni}
                    #return a dictionary with locus tag as index, and protein IDs as values
                    gene_prot_dict[locus_tag]=protein_ids
    else:
        print ('No UniProt accessions in:', genome_acc)
        
    return gene_prot_dict

def gene_centric_annotations(TF_instance):
    if not TF_instance.GO_term:
        return []
    go_annotations = []
    exp_regulations = get_exp_regulations(TF_instance)
    regulated_genomes = set(reg.gene.genome for reg in exp_regulations
                            if reg.expression_experimental_techniques)
    
    for genome in regulated_genomes:

        #get the mapping between genes in genome and uniprot accessions
        #this looks like:
        #{'CCNA_03751': {'refseq': 'YP_002519124.2', 'uniprot': u'A0A0H3CFU0'},
        #'CCNA_03752': {'refseq': 'YP_002519125.3', 'uniprot': u'A0A0H3CCD6'}}
        uni_refseq_map = uniprot_to_refseq(genome.genome_accession)
        
        if len(uni_refseq_map)>0:
            #get experimental regulations in genome [i.e. regulated genes]
            regulations_in_genome = exp_regulations.filter(gene__genome=genome)
    
            #for each experimentally-validated regulated gene instance
            for reg in regulations_in_genome:
                locus_tag = reg.gene.locus_tag   
                
                #if locus_tag has associated uniprot
                if locus_tag in uni_refseq_map:
                    #generate annotations
                    #for each evidence source for that regulation
                    for tech in reg.expression_experimental_techniques:
                        with_field='' #define WITH field, fill it up only if IPI evidence used
                        if tech.EO_term in get_ECO_IPI_terms():
                            with_field='RefSeq:'+reg.curation_site_instance.site_instance.genome.genome_accession
                        
                        #get uniprot accession
                        uniprot_id=uni_refseq_map[locus_tag]['uniprot']
    
                        print("Gene centric annotation:",uniprot_id)
                        
                        #make annotation
                        go_annotations.append(new_gene_centric_gpad_entry(
                            uniprot_id, TF_instance.GO_term.GO_term_id,
                            reg.ref_pmid, tech.EO_term, with_field))            
            
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
