"""This file contains titles and descriptions for curation steps."""

titles = {
    '0': "Publication selection",
    '1': "Genome and TF information",
    '2': "Experimental methods used in this paper",
    '3': "Reported sites",
    '4': "Exact site matches",
    '5': "Inexact site matches",
    '6': "Site Anntotation",
    '7': "Gene regulation (experimental support)",
    '8': "Curation information"
}
descriptions = {
    '0': "Please choose a publication to curate.",
    '1': """This step collects information on the transcription factor (TF), the
    specific strains reported in the manuscript and the NCBI GenBank sequences
    that reported sites and TF will be mapped onto.""",
    '2': """Select experimental techniques used to verify binding/expression of
    the sites reported in the curation. Provide a summary of the basic
    experimental procedure used to demonstrate binding/expression""",
    '3': '', # filled dynamically
    '4': """For each reported site, all exact matches in the chosen genome are
    listed. If a reported site does not have any exact matches, or the matched
    position/genes do not coincide with reported positions/gene, select the \"No
    valid match\" option. This will initiate a non-exact search.""",
    '5': """Inexact matches for sites without valid matches are listed here,
    sorted by affinity to the TF-binding motif.  If the matched position/genes
    do not coincide with reported positions/gene, select the \"No valid match\"
    option.""",
    '6': """Fill in the information regarding each site.""",
    '7': """Nearby genes are displayed for identified sites. Check all genes for
    which TF-site mediated regulation is reported in the manuscript. Skip this
    step if manuscript does not report gene expression.""",
    '8': """This step finalizes the curation. Fill all required fields."""
}
