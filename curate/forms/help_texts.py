"""This file contains all help texts that appear somewhere during paper
submission and curation processes. The purpose of keeping long help texts in a
separate file is to keep the code clean and make it easy to change those texts
later, if required."""

pubmed_publication_form = dict(
    pmid = "Paste the PubMed ID obtained from the NCBI website.",
    reported_TF = "Type the name of the transcription factor(s) reported in the manuscript.",
    reported_species = "Type the name of the species reported in the manuscript.",
    contains_promoter_data = "The paper provides experimental data on the structure and sequence of TF-regulated promoter.",
    contains_expression_data = "The paper provides experimental support for TF-mediated regulation of genes.",
    submission_notes = """Include any additional
details about the submission. For instance, you might indicate the approximate
number of sites reported, whether high-throughput techniques are used or any
other factor that might help prioritize curation.""",
)

non_pubmed_publication_form = dict(
    reported_TF = "Type the name of the transcription factor(s) reported in the manuscript.",
    reported_species = "Type the name of the species reported in the manuscript.",
    contains_promoter_data = "The paper provides experimental data on the structure and sequence of TF-regulated promoter.",
    contains_expression_data = "The paper provides experimental support for TF-mediated regulation of genes.",
    submission_notes = """Include any additional details about the
submission. For instance, you might indicate the approximate number of sites
reported, whether high-throughput techniques are used or any other factor that
might help prioritize curation.""",
)

# Curation help texts
publication_form = dict(
    pub = '',
    
    no_data= """Check this button if, after examining the paper, you find that
    the paper does not have data on binding sites. Checking this button will
    mark the paper as having no binding site data and set it to the 'curation
    complete' status. Also, the curation process will be ended as the paper has
    no data to be curated.""",
)

genome_form = dict(
    TF = """Select the transcription factor you are curating on from list. If not in list, please
    contact the master curator.""",
    
    TF_type = """If specified in the manuscript, select the quaternary structure for the transcription
    factor when binding to the sites reported in this curation.""",
    
    TF_function = """If specified in the manuscript, select the mode of operation for the TF on the sites
    reported in this curation.""",

    genome_accession = """Paste the NCBI GenBank genome accession number for the species closest to the
    reported species/strain.""",

    TF_species_same = '',

    site_species_same = '',

    TF_accession = """Paste the NCBI TF protein accession number for the species closest to the reported
    species/strain.""",

    TF_species = """Type the full name of the species/strain the TF belongs to as reported in the
    manuscript.""",

    site_species = """Type the full name of the species/strain in which the sites are reported in the
    manuscript.""",

    contains_promoter_data = "The paper provides experimental data on the structure and sequence of a TF-regulated promoter",

    contains_expression_data = "The paper provides experimental support for TF-mediated regulation of genes",
)

techniques_form = dict(
    techniques = """Select as many as apply to sites reported in this submission. Hover over any technique to see the
    description.""",

    experimental_process = """Write a concise, intuitive description of the experimental process to ascertain
    binding/induced expression""",

    external_db_type = """Select type of external database containing data (e.g. DNA-array data) reported in paper""",

    external_db_accession = """Type the accession number for external database referenced in paper.""",

    forms_complex = """Check if the manuscript reports characterization of the interaction of the TF with another
    protein""",

    complex_notes = """Provide brief description of the proteins involved in the complex and how it affects
    binding"""
)

site_entry_form = dict(
    is_motif_associated = """Check this option if sites are reported by the authors as being associated with
    a known TF-binding motif. Uncheck if sites are reported solely as DNA fragments shown to be bound by the TF,
    without reporting any specific binding sequences.""",

    sites = """Enter the list of sites in FASTA format or type the list of either site
    sequences or coordinates (one site per line). The sites can be entered in two
    major formats: sequenced-based (e.g. CTGTTGCACGT) or coordinate-based
    (e.g. 12312 12323). Optionally, quantitative data (q-val) can also be added
    to either format. All fields (i.e. site & q-val or coordinates & q-val) must
    be either space or tab separated.""",

    quantitative_data_format = """Please enter the data format for the quantitative values associated with sites""",
)

site_exact_match_form = dict()

site_soft_match_form = dict()

site_annotation_form = dict()

gene_regulation_form = dict()

curation_review_form = dict(
    revision_reasons = """Select, if needed, the reason why this curation may require revision.  See detailed list
    of reasons in the curation guide.""",

    confidence = """Check this if experimental techniques and results meet the
    standards specified in the curation guide""",
    
    NCBI_submission_ready = """A curation is ready for submission if: (a) the identified genome sequence
    matches the reported one or (b) identified and reported genomes match at the species level and at least
    90% of reported sites are located as exact matches.""",

    paper_complete = """Check this box if there are no more curations pending for this paper (additional
    sites, sites supported by different techniques, sites for other TFs, etc.""",

    notes = """Type in any additional notes on the curation process. For instance, if reported sites
    were left out for some reason, what prompted selection of a surrogate genome instead
    of another, general comments on the experimental process, etc.
    """,

    confirm = "Check to submit when you click \"next step\"",
)
