# Database

In CollecTF, all the collected data is stored in a relational database. Below, the
overall database structure is shown and tables are described.




-   `Curation` contains the data about the curation itself, such as the species of
    reported TF and sites, the overall experimental process, whether the TF is
    shown to interact with other protein/ligand that influences binding
    (i.e. does the TF form complex?). It also contains the meta-data about the
    curation, such as notes, if it requires revision and if it is validated by
    the master curator.
    
    Curation table has the following relations to other tables.
    
    -   `Curator`: The user who submitted this curation.
    
    -   `Publication`: The curated manuscript.
    
    -   `TF`: The transcription factor that is reported in the paper.
    
    -   `TF instances`: The specific instances of the TF reported in the paper. For
        example, if the TF is LexA and the paper is reporting binding sites in
        *E. coli*, the TF-instance could be the LexA protein in *E. coli* with the
        accession number NP<sub>418467</sub>.1. There is many-to-many relation between
        `Curation` and `TF instances` to be able to define TFs composed of several
        sub-units with different accession numbers (e.g. hetero-dimer such as IHF,
        composed of IHF alpha and beta).
    
    -   `Site instances`: The table that links the `Curation` to `Site
            instance` objects which are specific sequences on the genome, through an
        intermediate table, `Curation-Site-instance`.
    
    -   `ChIP info`: When the paper reports high-throughput data, the assay
        conditions and the high-throughput protocol is captured and linked to the
        `Curation`.
    
    -   `Quantitative data format`: If the manuscript reports quantitative values
        associated with each site, the format for that data is stored here. It can
        be ChIP-Seq fold enrichment, fluorescence ratio of ChIP-chip, etc.

-   `Curator` contains the curator type and user data. Curator type can be either
    `internal` or `external`. User data is stored using Django built-in
    user-management module and is one-to-one related to the `Curator` table.

-   `Publication` table keeps manuscript information such as PMID, authors,
    title, etc. It has also fields on whether it has promoter/expression data and
    reported TF/species that might be useful when the paper is being curated.

-   `Gene` table contains all genes in all genomes that is in the database. Each
    `gene` object contains position and strand information, locus tag and
    accession number.

-   `Genome` contains the accession number, organism information and key to the
    genome sequence that is stored in another table, due to performance
    issues. `Genome` table also has relation to the `taxonomy` which stores the
    taxonomy of organisms of which genomes are stored.

-   `Genome Sequence` contains only the genome sequence. In the initial
    implementation, it was a member of `Genome` table, however, it has been moved
    to a separate table, so that the whole genomic sequence is not fetched from
    database every time the `Genome` object is needed.

-   `Taxonomy` table stores the phylogeny in the database. It has the `rank`,
    which can be phylum, class, order, family, genus or species. Each taxonomy
    object has a relation to an object in the same table, to its parent. For
    instance, the object with the name *Escherichia coli* and rank *species* has
    link to its parent which is a `Taxonomy` object as well, with the name
    *Escherichia* and rank *genus*.

-   `TF` table contains the information about TF and a link to its family. For
    instance, the TF FNR is represented as a `TF` object with the name "FNR" and
    a link to a `TF Family` database entry with the name "FNR/CRP".

-   `TF Family` contains the name and description for the family.

-   `TF Instance` contains the information about specific TF proteins. Fields are
    *protein accession*, *name* and *description*. In the database design, the `TF
      instance` and `TF` tables are linked through `Curation` table, it would make
    more sense to link them directly, though.

-   `Site Instance` table stores the binding sites which are just specific
    sequences on the genome that have *start*, *end*, *strand* information and a
    link to the *genome* table. For getting the sequence of a binding site, the
    *sequence* field has been used, but due to data redundancy, it has been
    replaced with the method that extracts the genomic sequence from the genome,
    given the start and end positions. The *sequence* field is not used at all and
    it is kept just for sanity-check purposes.

-   `Curation-Site Instance` is the central part of the database. It links the
    `Curation` and `Site instance` tables with many-to-many relationship. In
    other words, a `Curation` might have multiple `Site instance` objects (i.e. a
    paper may report multiple binding sites) and a `Site instance` might have
    multiple `Curation` objects associated (i.e. the same binding site may have
    been reported by several papers). This table also contains additional
    information which are describe below.
    -   `site type` can be one of the following.
        -   *motif-associated* which is the exact binding site that the TF is shown
            to bind.
        -   *non-motif-associated* which is the region that the TF is shown to bind
            *somewhere* in it, but not specifically *where*.
        -   *variable-motif-associated* which the exact sequence that the TF is shown
            to bind, but can not be included in a gapless alignment with other binding
            sites. Binding sites that have variable spacing, inversion (direct-repeat
            in the motif and inverted-repeat in the specific binding site) fall under
            this category.
    
    -   `annotated sequence` is the *actual* sequence that is reported in the
        paper. The whole curation process can be seen as the mapping of binding
        sites reported in the paper to the genome in the CollecTF. In most cases,
        the genome that has been used in the experiment is available in the NCBI
        RefSeq database, so it can be used to find the binding sites in the
        genome. However, if the strain that has been used is not fully sequenced
        and annotated, the curator selects the closest genome in the RefSeq in
        terms of phylogeny. In this case, the exact reported site may not be found
        in the genome (because the genome belongs to some other
        species!). If so, the `annotated sequence` stores the reported binding site
        where the `Site instance` object itself has the data for the surrogate
        binding site.
    
    -   `TF type` is the structure of the TF (e.g. monomer, dimer, tetramer
        ,etc.). Initially, this field was part of `Curation` table. However,
        although unlikely, it is possible for a specific TF to have different
        structures for different binding sites. Initial design would require
        multiple curations in such a case. It has been moved to this table so that
        sites that TF binds with different structures can be submitted within the
        same curation.
    
    -   `TF function` can have values *activator*, *repressor*, *dual* or *NA (not
        specified)*. Like `TF type` it has been part of `Curation` table. It has
        been moved into this table to allow a single `Curation` to have associated
        site instances that the TF is acting as a repressor and site instances that
        the TF is acting as an activator.
    
    -   `experimental techniques` is a many-to-many relation to `Experimental
            technique` table. It stores the list of experimental technique that have
        been used to identify the binding site.
    
    -   `regulates` is a many-to-many relation to the `Gene` table through
        `Regulation` table which stores the information of which genes are
        inferred/verified to be regulated by the TF by binding to the site.
    
    -   `quantitative value` is an optional field. If the paper reports sites with
        associated numerical values (e.g. ChIP peak intensities), this is the field
        that the value is stored
    
    -   `is-obsolete` and `why-obsolete` are NCBI submission-related fields. If a
        binding site has been identified as incorrect due to some reason and if it
        has been already submitted to the NCBI RefSeq database, it can not be simply
        removed from the database as it will result to dead DBXRef in the RefSeq
        page. Instead, it will be marked as obsolete with the note explaining why it
        has become obsolete. When it is requested through the DBXRef link, an
        appropriate message will be presented.

-   `Regulation` stores the TF regulation of genes. It links two tables:
    `Curation-Site instance` and `Gene`, to capture which TF regulates which gene
    by binding which site on the genome. In addition, the evidence type is
    stored. If there is experimental evidence of regulation in the paper, saying
    that TF *x* up/down regulates gene *y*, the evidence type is
    `exp-verified`. Otherwise, all genes in the operon of which the site is
    upstream of are labeled as `inferred`. See operon prediction.

-   `Not Annotated Site Instance` has two fields: `sequence` and key to the
    `Curation` table. For some curations, the binding site that is reported in
    the paper can not be matched to any sequence in the genome for some reason
    (e.g. the reported genome is different from the genome used in the
    paper). Therefore, it is stored as `Not Annotated Site Instance`.

-   `Experimental Technique` is the table to store experimental techniques that
    are used to identify TF-binding or regulation. Each binding site is
    identified experimentally using a set of techniques. This table has fields
    for the technique name, description and type which specifies if it is used to
    show binding or expression or if it is an in-silico technique.

-   `Experimental Technique Category` can be considered as labels that define
    experimental techniques. It is different from the `type` field in
    `Experimental Technique` table as a technique may have many categories at the
    same time. For instance, the technique *ChIP-PCR* has three categories:
    *PCR-based techniques*, *Immunoprecipitation methods* and *in-vivo
    techniques*.

-   `ChIP Info` stores the information about curations that have high-throughput
    data submission. If the paper reports binding sites using high-throughput
    methods, the assay conditions and method notes are stored in this table which
    is linked to the `Curation` table.

-   `External Database`. Sometimes, authors choose to upload additional data
    (e.g. DNA-array or gene-expression data) to a database. The information about
    those databases are stored in this table. For each external database, the
    name, description and URL format is stored.

-   `Curation-External Database` links the `Curation` table to `External Database`
    table, if there is any data in the paper that is available in an external
    database. If so, the fields `accession number` and key to `External Database`
    enable to access to that piece of data.

-   `NCBI Submission`. The curated data in CollecTF is integrated into NCBI RefSeq
    database through periodic genome-specific submissions. The internal table
    keeps track of all the data that has been submitted to the NCBI RefSeq so
    far. For each `Curation-Site instance` that has been submitted, there is an
    entry in this table having *submission time* and the *genome accession number*
    that the data is submitted to. The reason for explicit *genome accession
    number*, which can be inferred through the key to `Curation-Site instance` is
    to be able to handle major genome updates in NCBI part (e.g. `NC_000913.2` to
    `NC_000913.3`). Although we don't have a specific strategy how to modify all
    CollecTF records, the explicit *genome accession number* field allows us to
    keep track of which genome that the NCBI submission is made for, even if we
    decide to update `Genome` relations of all `Site instance` objects from
    `NC_000913.2` to `NC_000913.3`.
