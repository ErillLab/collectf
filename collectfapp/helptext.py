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
    'contains_promoter_data': "The paper provides experimental data on the structure and sequence of a TF-regulated promoter",

    'contains_expression_data': "The paper provides experimental support for TF-mediated regulation of genes",

    'techniques':  """Select as many as apply to sites reported in this submission. Hover over any technique to see the
    description.""",

    'experimental_process': """Write a concise, intuitive description of the experimental process to ascertain
    binding/induced expression""",

    'external_db_type': """Select type of external database containing data (e.g. DNA-array data) reported in paper""",

    'external_db_accession': """Type the accession number for external database referenced in paper.""",

    'forms_complex': """Check if the manuscript reports characterization of the interaction of the TF with another
    protein""",

    'complex_notes': """Provide brief description of the proteins involved in the complex and how it affects
    binding"""

}

site_report_form = {
    'is_motif_associated': """Check this option if sites are reported by the authors as being associated with
    a known TF-binding motif. Uncheck if sites are reported solely as DNA fragments shown to be bound by the TF,
    without reporting any specific binding sequences.""",

    'is_coordinate': """Check this option to enter reported sites as coordinates instead of sequences.
    Coordinates should be {start end}. Coordinate order determines strand. One site per line.""",

    'is_chip_data': """Check if reported sites come from ChIP experiments (e.g. ChIP-chip, ChIP-seq).
    A few ChIP-specific fields will be populated.""",

    'has_quantitative_data': """Check to report quantitative data associated with individual site instances (e.g.
    binding affinity). Quantitative data for each site may be entered at the end of the line, separated by tab/space.
    For sequence entry: {sequence value}. For coordinate entry: {start end value}.""",

    'peak_calling_method': """Describe briefly the nature and range of quantitative data (e.g. intensity readout from
    EMSA. Range: 2.3 - 5.7). For ChIP experiments, simply state: 'ChIP peak intensity' and list the peak calling method
    used.""",

    'assay_conditions': """Describe the conditions of the ChIP experiment that capture the specifics of the in-vivo 
    setting for cross-linking. Were cells at exponetial-phase? Was the system induced? How were cells grown?""",

    'method_notes': """Describe (use copy-paste if appropriate) the ChIP protocol. What antibodies were used? What
    chip/sequencer and using what parameters? Etc.""",

    'sites': """Type either site sequences or coordinates. One site per line.
    FASTA format is supported for sequence entry mode.""",

    'chip_data_extra_field': """Quantitative data can be automatically associated to motif-associated sites. Paste 
    corresponding ChIP peak data here {peak_start_pos peak_end_pos intensity_value}. The system will search entered sites
    in peaks and associate peak intensities. Mapped peak intensity values will be displayed for review before curation submission.""",

}

site_exact_match_form = {
    
}

site_soft_match_form = {
    
}

site_regulation_form = {
    
}

curation_review_form = {
    'revision_reasons': """Select, if needed, the reason why this curation may require revision.  See detailed list
    of reasons in the curation guide.""",

    'confidence': """Check this if experimental techniques and results meet the
    standards specified in the curation guide""",
    
    'NCBI_submission_ready': """A curation is ready for submission if: (a) the identified genome sequence
    matches the reported one or (b) identified and reported genomes match at the species level and at least
    90% of reported sites are located as exact matches.""",

    'paper_complete': """Check this box if there are no more curations pending for this paper (additional
    sites, sites supported by different techniques, sites for other TFs, etc.""",

    'notes': """Type in any additional notes on the curation process. For instance, if reported sites
    were left out for some reason, what prompted selection of a surrogate genome instead
    of another, general comments on the experimental process, etc.
    """,

    'confirm': "Check to submit when you click \"next step\""
}
