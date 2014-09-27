
# Browse and Search CollecTF

CollecTF allows browse the database by TF, species and experimental techniques
as well as fully customizable search where any subset of TFs, species can be
searched for trasncription factor binding sites that are identified by any of
the selected experimental techniques.

# Browse CollecTF

CollecTF can be browsed at three different levels

-   NCBI taxonomy
-   TF or TF family
-   Experimental technique

## Browse by NCBI taxonomy

CollecTF has species information for each binding site record, as well as all
levels of taxonomy, up to the bacteria domain. Browse page enables user to
select any level of taxonomy to see the curated binding sites. For any selected
level in the taxonomy, the related information is fetched from Wikipedia, the
database is searched by the selected criteria and all results are listed,
grouped by TF and species. The list has links to individual reports, each
containing the set of binding sites for a particular TF and species. The result
page also has the link to the combined report which enables the user to see all
individual reports in one page.

## Browse by TF and TF family

For a selected TF or TF-family, all binding site records are retrieved from the
CollecTF database. On top of the result page, the description of the TF or
TF-family is displayed, retrieved from the database. Below that, the list of
reports is displayed which contains links to individual reports where the user
can see all binding sites, the motif, etc.

## Browse by experimental technique

The user can also search binding sites that are identified using a particular
experimental technique or a technique type.

# Search CollecTF

Search page provides full customization on the filter that is used to search
and browse the database. The search consists of three steps:
1.  Select transcription factor(s): any number of TFs or TF families can be
    selected. The search box on the list provides a quick way to find the
    TF/family of interest.
2.  Select species: any number of species (or classes, orders, etc.) can be
    selected. Similar to the previous step, the search box can be used to find
    the taxonomy unit easily.
3.  Select experimental techniques: Any combination of experimental techniques
    can be selected. There are three selection sets are provides which can be
    linked using Boolean operators.

Based on the selected options, the database is searched for binding sites
satisfying the selected criteria and the results are listed, grouped by TF and
species. Using the result list, the user can access the collection of binding
sites of a particular TF and species.

# How are report pages generated?

Based on selected search/browse criteria, all collection of binding sites
(`Curation-SiteInstance` objects) are filtered. The filtered collection may
contain binding sites from multiple TFs and species. After grouping them by TF
and species, the TF-species-specific binding site collection is processed to
generate the report page. This section describes the report generation process.

The method `make_reports` in `/collectf/browse/motif_report.py` is called with
all filtered `Curation-SiteInstance` objects. They are grouped by TF and
species. Furthermore, each TF-species-specific group is split into two:
motif-associated and non-motif-associated `Curation-SiteInstance` objects. For
each group (pair of motif/non-motif associated objects), the `MotifReport`
object is created where the meta-sites are computed.

## What is meta-site?

Experimental evidence of binding for a given genomic position may be distributed
across multiple reported sites of either type (i.e. motif-associated and
non-motif-associated). Wherever a motif-associated site has been defined,
CollecTF dynamically combines multiple sources of evidence by arbitrarily
defining a leader site and using two simple pair-wise propagation
rules. Evidence from two motif associated sites is combined if the overlap
between sites is >75% of the combined site length. Evidence from non-motif
associated sites is integrated into a leader site if they fully overlap any of
the combined motif-associated sites. See `/collectf/base/metasite.py` for
details.

## Site alignment

The generated report also has the alignment of binding sites. Since the sites
for a given TF and species might be reported in different lengths, a naive
sequence alignment would not be helpful. To produce alignments of variable
length sequences, [LASAGNA](http://www.ncbi.nlm.nih.gov/pubmed/23522376) algorithm is used.
