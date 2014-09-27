# Curation

# Introduction

This section of the CollecTF documentation covers the curation process, specifically
-   the content of forms that are filled and submitted by the curator and
-   how they are processed in the background.

# Data Submission

## New Publication Submission

Each curation in the CollecTF is associated with a publication. The first step
on the process of curating the data in the paper to the CollecTF database is to
add the paper to the collection to be curated later.

The link to add paper is available [here](http://collectf.umbc.edu/curate/add). To add the paper, the user must be
logged in. When the page is requested, the user is asked to provide the
publication information via submission of a form having the following fields.
-   **PMID:** The PubMed ID for the paper. It can be obtained from the NCBI
    website. It is used to retrieve information about the paper
    (e.g. title, authors, the journal and the year it has been
    published) from NCBI database (through [NCBI Entrez](http://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/))
-   **Reported TF(s) and reported species:** The transcription factor and the
    species reported in the manuscript.
-   **Promoter/Expression information:** Two check boxes about whether the
    manuscript contains promoter and expression information.
-   **Submission notes:** Anything that may be useful for the curation.

When the form is submitted, it is validated, the paper information is retrieved
from NCBI Publication database and a new form containing paper-related
information (e.g. title, authors, journal, etc.) is presented to the user for
confirmation. Besides confirmation, the user can
-   mark himself/herself as the to-be-curator of the paper, and
-   check the paper as having no data. If the paper has no data on binding sites,
    this option is selected and the paper is marked as 'complete' and 'having no
    binding site data'.

For the implementation details, see `curate/add_publication.py` and
`curate/forms/add_publication_form.py`.

## Curation Submission

Curation submission is the most important component of the whole database and
web interface. All the binding site data, that the CollecTF is all about, is
submitted through the curation submission form. The presentation and back-end
processing of the curation submission form is split into many components which
are explained in this section.

The submission pipeline for conventional data and high-throughput data
(e.g. reported binding sites together with peaks) are slightly different. In the
next section, submission pipeline is described for conventional data
submission. In the later section, the differences for high-throughput data
submission are described.

### Submission Pipeline

The curation submission process contains of a few steps where each step consists
of a form that is submitted along with the required information. Until last
step, all of the entered data are stored in session structure. The session
structure is also useful to process forms that are interdependent to other
steps.

The curation submission page is available [here](http://collectf.umbc.edu/curate/curation/), only for registered users. Here
is how the pipeline works. When the user requests this page, the first step of
the pipeline is presented to the user. The user enters the required data. The
form is checked if it contains everything required. If not, the same form is
presented again, with the messages pointing the problem. If the form is valid,
the entered data are processed, stored as session data, if necessary. The next
step is identified and presented to the user. When the user submits the form in
the last step, all the data entered are checked again for consistency. If
everything seems valid, the data is saved into the database and the process is
completed.

1.  Submission Pipeline Implementation Details

    All the source code related to the curation submission is under `curate/`
    directory. The submission pipeline is implemented using `form-wizard`
    application in `Django`.
    
    From the official `Django` documentation, here is how the `form-wizard` works.
    
    -   The user visits the first page of the wizard, fills in the form and submits it.
    -   The server validates the data. If it’s invalid, the form is displayed again,
        with error messages. If it’s valid, the server saves the current state of the
        wizard in the backend and redirects to the next step.
    -   Step 1 and 2 repeat, for every subsequent form in the wizard.
    -   Once the user has submitted all the forms and all the data has been validated,
        the wizard processes the data - saving it to the database, sending an email,
        or whatever the application needs to do.
    
    More documentation on `Django` `form-wizard` is available [here](https://docs.djangoproject.com/en/1.6/ref/contrib/formtools/form-wizard/).
    
    The entry point for the submission pipeline is the `curation` function in the
    `add_curation.py` file. This method also contains the list and order of the
    forms of the pipeline. The `curate/forms/add_curation_form.py` contains the
    definition of each form (i.e. step of the pipeline). Each form is implemented as
    a class which contains all of the fields of the form, as well as methods to
    validate them.
    
    Before a form is presented, its `get_form` function is called which constructs
    the form for the given step. When the form is submitted, `clean` functions
    defined in the form class definition are called. If the form is valid, the
    `process_step` function of the corresponding step is called, to process entered
    data and store it in the session before proceeding to the next step.

2.  Publication Selection

    The first step of the submission pipeline is the selection of the paper to be
    curated. The papers that are already assigned to the user are retrieved from the
    database and listed for selection. When the curation is submitted, the paper is
    linked to the curation as a link between `Publication` and `Curation` tables.
    
    This step also contains a check box which says "This paper contains no
    data.". This option is checked if, after examining the paper, the curator finds
    that the paper does not have data on binding sites. Checking this check box
    marks the paper as having no binding site data and set it to
    `curation-complete` status. Also, the submission pipeline is terminated as the
    paper has no data to be curated.

3.  Genome and TF Information

    After the selection of the paper to be curated, the curator is asked to enter
    the genome and TF information. Specifically, the curator is asked to fill the
    following fields in the form.
    
    -   `TF`: This is a drop-down list containing all the transcription factors in the
        database. The curator selects the TF that the paper reports binding sites
        of. If a specific TF is not on the list, the curator can go to
        <http://collectf.umbc.edu/curate/add_TF/> to add the TF to the list.
    
    -   `Genome NCBI accession number`: The genome accession number of the genome is
        entered here. This field is a text-box which requires the genome accession
        number with the version number (e.g. `NC_000913.2`, not
        `NC_000913`). Auto-complete feature on this field tries to complete the genome
        accession number based on genomes that are already in the database. For
        instance, as the curator starts typing `NC_000913.2`, all genomes having an
        accession number that starts with `NC_..` are displayed. The auto-complete
        feature works for species names as well. For instance, if the curator starts
        typing `Escherichia coli`, only the accession numbers of genomes *Escherichia
        coli* strains are appeared.
        
        Since an organism may have more than one chromosome and some number of
        plasmids, a paper may report binding sites for more than one
        chromosome/plasmid sequences where each of them has their own genome accession
        number. To handle such cases, extra genome accesssion fields are used to enter
        genome accession numbers if needed. Extra genome accession fields are hidden
        by default can be toggled if necessary.
        
        *Implementation note*: If the genome of the entered accession number is
        already in the database, it is fetched from the database. If not, the genome
        accession number is used to connect to NCBI genome database to download the
        `genbank` file for that genome. The file is then parsed, it is stored in the
        `Genome` and `GenomeSequence` tables of the database. The genes are also
        parsed from the `genbank` file, all genes are stored in the database as well,
        as entries in the `Gene` table.
        
        The number of `Genome accession number` fields is defined in the CollecTF
        configuration file, sitting in `collectf` directory.
    
    -   `TF accession number`: The NCBI TF protein accession number is entered
        here. Like `genome accession number` field, this field has auto-complete
        feature which tries to complete the accession number based on TF proteins in
        the database. The auto-complete feature can also be used with the species name
        (e.g. by typing *Escherichia coli*, all TF proteins for *E. coli* in the
        database are displayed). The accession number doesn't need to have the version
        number (e.g. `NP_418467`, not `NP_418467.1`). If it has, it is removed
        automatically before the database look-up and retrieval from NCBI Protein
        database if it is not already in the CollecTF database.
        
        A TF may be composed of several sub-units where each sub-unit has its own
        accession number (e.g. IHF composed of IHF-alpha and IHF-beta). To handle
        such cases, more than one `TF accession number` field can be used. Extra
        fields are hidden by default, but can be toggled if needed.
        
        *Implementation note*: If the TF of the entered accession number is already
        in the database, it is fetched from it. Otherwise, the TF accession number is
        used to connect to NCBI protein database to download the record which is used
        to create `TF Instance` object in the database.
        
        The number of `TF accession number` fields is defined in the CollecTF
        configuration file, sitting in `collectf` directory.
    
    -   `Site species` is the text field that should contain the full name of the
        species/strain in which the sites are reported in the manuscript. Since the
        genome of which the accession number is entered can be different from what is
        reported in the manuscript (i.e. the reported genome may not be available in
        the RefSeq database), this additional field is used to get the full name of
        the species/strain that is reported in the paper.
        
        If the reported strain is the same with the genome, the check box saying
        "This is the exact same strain reported in the manuscript for the sites" is
        checked which makes the `Site species` text box disabled.
    
    -   `TF species` is the text field that should contain the full name of the
        species/strain the TF belongs to as reported in the manuscript. Like `Site
          species` field, it is necessary to have the strain/species of the TF used in
        the paper since the TF protein entered in the `TF accession number` field can
        be different from what is reported in the paper.
        
        If the accession number of the TF that the paper reports is used, the check
        box saying "This is the exact same strain as reported in the manuscript for
        the TF" is checked which makes the `TF species` text box disabled.
    
    -   `Contains promoter/expression data` fields are used to indicate whether the
        paper has promoter and expression data. These two fields are also used in
        "publication submission" form. The `Publication` object is modified based on
        these two check boxes. Also, if the `contains expression data` check box is
        not checked, the *Gene Regulation* step becomes read-only.

4.  Experimental Methods

    In this step, *all* experimental techniques used to verify binding/expression
    of the sites are selected. In addition, the curator is asked to provide a brief
    summary of the experimental procedure.
    
    -   `Techniques` is the list of techniques available in the database.
    -   `Experimental process` is the text field that should contain the summary of
        the experimental process.
    -   `Additional information` section contains two check-boxes.
        -   *The manuscript reports high-throughput data from an external database.* If
            checked, extra fields appear where the curator can provide external
            database accession numbers if applicable (e.g. GEO accession number).
        -   *The manuscript reports that TF forms complex with other proteins for
            binding with reported sites.* If checked, a text-box appear where the
            curator can provide detail on the TF forming complex with other molecules
            in order to bind.

5.  Site Entry

    In this step, the curator enters the binding sites reported in the paper.
    
    -   `Site type` TF-binding sites can be defined at different levels. By
        definition, a TF-binding site is simply a (relatively short) stretch of DNA to
        which a transcription factor is shown to bind (e.g. a ChIP-Seq peak). Many TFs
        target known specific sequence patterns in the DNA. Some of these patterns are
        complex and require gapped alignment alignment (e.g. because of variable
        spacing) or more complex procedures in order to be defined. Other patterns are
        simpler and can be represented by a gapless alignment of sites (known as a
        motif), providing a much more concise definition of TF-binding site. In
        CollecTF we refer to these site types as *motif-associated* [for gapless
        alignment], *variable motif-associated* [for complex patterns] and *non-motif
        associated* [for unknown or absent patterns; just evidence of binding].
    
    -   `Sites` is text-box for sites. The sites can be entered in FASTA
        format or list format. In the list format, the sites can be entered in two
        different ways: sequence-based (e.q. `CACCACACGATCGATC`) or coordinate-based
        (e.g 12312 12323). Optionally, quantitative data (`qval`) can also be added
        to either format. All fields must be either space or tab separated.
    
    -   `Quantitative data format`. If the `sites` field contains quantitative values
        (e.g. peak intensities, estimated K<sub>d</sub>), the format for these values is put
        into this field.

6.  Exact Site Matches

    For the sites entered in the previous step, there may be multiple instances in
    a genome. After submission, sites submitted as sequences are searched in the
    genome, specified in `Genome and TF information` step. Exact matches to
    submitted sites are reported back specifying their location in the genome and
    nearby genes. For a given site, if there are more than one matching, the
    curator is asked to select one. For sites that have no exact match, the only
    option will be "No valid match". Also on the top-right of the page, the
    sequence logo is displayed, built from the sites that are entered in the
    previous step.

7.  Inexact (soft) Site Matches

    In some cases, especially if using a sequence that is not an exact match to the
    reported strain, some sites may not be found using an exact search. In this
    case, CollecTF uses the available evidence to construct scoring matrix and
    search the genome for slightly inexact matches (up to 2 mismatches away from the
    reported site). 
    
    Note that, only the sites that are can not be matched in `Exact site match`
    step, are searched in this step.
    
    *Implementation note*: For a given site, the inexact matches are sorted by PSSM
    score where the matrix is built from the sites that are matched in the `Exact
    site matches` step. To see how exact and inexact matches are performed see
    `/collectf/curate/site_entry.py`. <span class="timestamp-wrapper"><span class="timestamp">&lt;2014-05-16 Fri&gt;</span></span> It seems that motif is not
    used to score the inexact matches. To be fixed.

8.  Site Annotation

    During this step, specific experimental techniques are matched to individual
    sites already identified in the reference genome. The quaternary structure of
    the TF when interacting with sites (e.g. dimer), as well as the regulatory mode
    on TF-binding at each site (e.g. repressor), if known, can also be entered
    independently for each site. In addition, if quantitative data for sites have
    been manually entered, they can be validated here.
    
    The user can select multiple sites using the mouse in combination with `Shift`
    key or through `Select/Unselect all` link to easily assign attributes to
    several sites at once.

9.  Gene Regulation (Experimental Support)

    If the paper is marked as having regulation information, the user is asked, for
    each reported site, which genes are shown to be regulated by the TF.

10. Curation Review and Submission

    -   `Revision required` If the curation is required to be revised for some reason
        (e.g. no comparable genome in NCBI, matching genome still in progress, etc.)
        the appropriate reason is selected.
    -   `I am confident of the results reported in this manuscript` If the
        experimental techniques and results meet the standards specified in the
        curation guide, this check-box is selected by the curator.
    -   `Curation is ready to submit to NCBI` If (a) the identified genome sequence
        matches the reported one, or (b) identified and reported genomes match at the
        species level and at least 90% of the reported species are located as exact
        matches.
    -   `Curation for this paper is complete.` If there are no more curations pending
        for the paper, this check-box is checked.
    -   `Notes` field is for any additional notes on the curation process.
    -   `I want to submit this curation` field must be checked before submitting the
        curation.

## Curation Submission for High-throughput Data

Curators can check ChIP and other high-throughput methodologies in the regular
submission mode, but if they are submitting data that is primarily based on
high-throughput binding assays (e.g. ChIP-seq, genomic-SELEX, etc.) they are
then encouraged to use the high-throughput submission portal.

The only difference over regular submission pipeline is the `site entry`
step. In this step, in addition to the `sites` box that sites can be entered,
curators are also able to enter peak data. This is most likely to be in
coordinate mode, but sequence entry is also accepted and detected
automatically. The sites in `sites` box are saved as site type selected
(i.e. motif associated, variable motif associated or non-motif associated) and
the peaks are saved into database as non-motif associated sites. If
quantitative values are provided for peaks and not for sites, the system
searches peaks for the sites and associate the retrieved peak quantitative
values to the sites.

The system also asks for assay conditions and ChIP methodology notes.

## Adding New TFs and TF Families

Many TFs and their families are already in CollecTF database. However, in case
of a paper reporting for a TF that is not in the database, the user is able to
add new TFs and its family, if necessary.

## Adding New Experimental Techniques

The curator can also add a new experimental technique, if it doesn't exist in
CollecTF. In addition to the technique name and description, technique type
(i.e. detection of binding, assessment of expression or in-silico prediction)
and a set of categories must be selected by the user.
