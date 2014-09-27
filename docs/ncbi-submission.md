# NCBI submission

CollecTF entries are periodically submitted to [NCBI](http://www.ncbi.nlm.nih.gov/) for integration into [RefSeq](http://www.ncbi.nlm.nih.gov/refseq/)
complete genome records as link-out features. For example, let's take a look at
[Vibrio cholerae O1 biovar El Tor str. N16961 chromosome I, complete
sequence](http://www.ncbi.nlm.nih.gov/nuccore/15640032). CollecTF entries are integrated as `protein_bind` features which
contain
-   location of the binding site (e.g. `74751..74771`)
-   `/experiment` field for the list of experimental techniques used to identify
    the binding site, with a link to the paper reporting it
    (e.g. `/experiment = ChIP-Seq [PMID:21750152 ]`)
-   `/note` for the feature description (e.g. `/note = Transcription factor
      binding site for NP_231738`)
-   `/bound_moiety`, the TF that binds
-   `/db_xref` link to the CollecTF page for the binding
    site. (e.g. `/db_xref = CollecTF:EXPSITE_00000230` pointing [here](http://collectf.umbc.edu/EXPSITE_00000230))

## `tbl` export

NCBI annotations are submitted in a form of NCBI feature table (`tbl` file),
see description at <http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2#tbl>. The `tbl`
file can be generated [here](http://collectf.umbc.edu/ncbi/export_tbl_view) only by a super user. When the form is submitted
with the genome accession number, the server returns a zip file containing the
files required for NCBI submission.

When the form is submitted, all `Curation-SiteInstance` objects for selected
genome are filtered. They are grouped by TF instances and meta-sites are built
from them. For each meta-site, if any of its sites has been submitted before, it
is skipped. Otherwise, its leader site is marked for submission and
corresponding `NCBI submission` object is created.

## Generating .sqn file

The NCBI submission process is still not fully automated. The command line
utility [tbl2asn](http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) has to be run by the superuser to generate the `sqn` file which
is uploaded to the NCBI ftp server. The downloaded archive file needs to be
extracted, the `tbl2asn` must be copied into that directory and run in that
directory by typing

    ./tbl2asn -p . -x .asn -Vvbr

which generates `.sqn` file for submission to Genbank.

## Database

CollecTF keeps track of what has been sent to NCBI RefSeq to avoid redundant
submissions. For each binding site annotation sent to NCBI RefSeq, an `NCBI
Submission` object containing
-   the submission time
-   the genome accession number (and the version number) that the curated data is
    sent for.
-   the `Curation-SiteInstance` object.

The "export tbl" page has a checkbox that should be used when the export is for
test purposes only. By checking the checkbox, no database change will be made.
