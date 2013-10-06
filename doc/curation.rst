Curation Process
================

Publication submission
----------------------
Publication submission is fairly simple. User puts the PMID for the paper. After
confirmation page, if the paper is not already in the database, it is added. The
paper can also be assigned to the curator.

Curation submission
-------------------

Curation submission consists of several steps. For each step, a form is displayed,
the curator fills all fields properly, submits the form and proceeds to the next
step.

- Publication selection :: All papers assigned to the curator are listed. The curator
     selects the paper that he/she is curating and moves to the next step.

- Accession numbers and general information :: This step contains TF and genome
     information reported in the paper.

- Methodology :: Experimental techniques used in the paper are submitted here.

- Site entry :: The reported sites are entered in this step.

- Exact site search :: Sites that are entered in the previous step are searched in
     the genome sequence. For a given site, there may be more than one exact
     matches. For each site that has at least one exact match are presented in this
     step for curator to match reported sites and matched instances in the genome.

- Soft site search :: For some papers, the reported genome is not available in RefSeq
     database, so the curator picks the closest genome from NCBI RefSeq database. In
     that case, for some reported sites, the exact match may not be present in the
     surrogate genome.
   
     For all reported sites that are not matched to any site instance in the previous
     step (exact site search), a soft search is performed to find site instances that
     are not exactly same but similar to the reported site (upto 2 mismatches).

Operon Prediction
-----------------
When the curator enters list of sites reported in the paper, all sites are searched
in the genome and all matches are returned. For each site-match in the genome, in
addition to the genomic position of the match, nearby genes are identified so that
they can be marked as regulated by TF (if any experimental evidence in the paper).

For each identified site instance, operons in both directions (if any) are
returned. The implementation for operon prediction is given below.

External Submission
-------------------
Curation submission form for external users is same as default form except few
changes in the last step (curation review step). For external users, the field
``revision_reasons`` is hidden and set to ``external_submission`` and the field
``NCBI_submission_ready`` is hidden and set to false.

Currently, user registration is disabled and is handled via feedback form (user requests
an account and we open a new account through admin page).


