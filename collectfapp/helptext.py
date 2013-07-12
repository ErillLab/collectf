"""This file contains all help texts appear somewhere during curation process. The
purpose of keeping long help text in separate file is to keep Python code file clean."""

publication_form = {
    'pub': '',
    
    'no_data': """Check this button if, after examining the paper, you find that the paper does not
    have data on binding sites. Checking this button will mark the paper as having no
    binding site data and set it to the 'curation complete' status. Also, the
    curation process will be ended as the paper has no data to be curated.""",
}

genome_form = {
    'TF': """Select the transcription factor you are curating on from list. If not in list, please
    contact the master curator.""",
    
    'TF_type': """If specified in the manuscript, select the quaternary structure for the transcription
    factor when binding to the sites reported in this curation.""",
    
    'TF_function': """If specified in the manuscript, select the mode of operation for the TF on the sites
    reported in this curation.""",

    'genome_accession': """Paste the NCBI GenBank genome accession number for the species closest to the
    reported species/strain.""",

    'TF_species_same': '',

    'site_species_same': '',

    'TF_accession': """Paste the NCBI TF protein accession number for the species closest to the reported
    species/strain.""",

    'TF_species': """Type the full name of the species/strain the TF belongs to as reported in the
    manuscript.""",

    'site_species': """Type the full name of the species/strain in which the sites are reported in the
    manuscript.""",    
}

techniques_form = {
    'contains_promoter_data': "The paper provides experimental data on the structure and sequence of TF-regulated promoter",

    'contains_expression_data': "The paper provides experimental support for TF-mediated regulation of genes",

    'techniques':  """Select as many as apply to reported sites. Hover over any technique to see the
    description.""",

    'experimental_process': """Write a concise, intuitive description of the experimental process to ascertain
    binding/induced expression""",

    'external_db_type': '',

    'external_db_accession': '',

    'forms_complex': '',

    'complex_notes': """Provide brief description of the proteins involved in the complex and how it affects
    binding"""

}

site_report_form = {
    'is_motif_associated': """If reported binding sites are 'true' binding sites, check this option. If a
    reported site is much longer than actual binding site sequence and contains the true binding
    site somewhere in its sequence, it is called a non-motif-associated site. In
    that case, this option must be un-checked.""",

    'is_coordinate': """If the paper reports coordinates instead of actual sequences, check this
    option. Enter coordinates for all sequences below (one for each line). Each line should have
    start and end positions, separated by space or tab characters. If the reported sites have quantitative data,
    they should be given as third item for each site instance (i.e. {start end quantitative_value}).""",

    'is_chip_data': """If the paper reports sequences/coordinates identified using ChIP methods, check this
    box. A few extra fields will be populated.""",

    'has_quantitative_data': """If the paper has some quantitative values associated with each site, please check
    this option. For example, if the used technique is ChIP, these values are peak
    intensities. Also, for some other tecnhiques, for each site, a quantitative value
    may be reported. They should be given as the last item for each site instance
    (i.e. If sequences are reported, the format would be {sequence value}, and if
    the coordinates are reported, the format would be {start_pos end_pos value}.""",

    'peak_calling_method': '',

    'assay_conditions': '',

    'method_notes': '',

    'sites': """Type either site sequences or coordinates. Site instances must be given in separate lines.
    If only sequence input is given, FASTA format is supported.""",

    'chip_data_extra_field': """If the used technique is ChIP and reported sites are motif-associated,
    paste ChIP peaks to this field with intensity values (the format is {peak_start_pos peak_end_pos intensity_value}).
    The identified motif-associated sites (given in 'Sites' field) will be searched in list of peaks (given in this field),
    and peak intensities will be associated with motif-associated sites. The mapped peak intensity values will be displayed
    before curation submission for review.""",

}

site_exact_match_form = {
    
}

site_soft_match_form = {
    
}

site_regulation_form = {
    
}

curation_review_form = {
    'revision_reasons': """Select, if needed, the reason why this curation requires revision.  See detailed list
    of reasons in the curation guide.""",

    'confidence': """Check this if experimental techniques and results meet the
    standards specified in the curation guide""",
    
    'NCBI_submission_ready': """A curation is ready for submission if: (a) the identified genome sequence
    matches the reported one or (b) identified and reported genomes match at the
    species level and at least 90% of reported sites are located as exact
    matches.""",

    'paper_complete': """Check this box if there are no more curations pending for this paper (additional
    sites, different techniques, other TF, etc.""",

    'notes': """Type in any additional notes on the curation process. For instance, if reported sites
    were left out for some reason, what prompted selection of a surrogate genome instead
    of another, general comments on the experimental process, etc.
    """,

    'confirm': "Check to submit when you click \"next step\""
}
