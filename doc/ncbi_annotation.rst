NCBI RefSeq Annotation
======================

Meta-site Rule
--------------
Experimental evidence for a given binding site may be reported in more than one
paper. In addition to that, in some cases, the coordinates for the same binding
location may be slightly different (a few bp). To avoid listing same binding
locations multiple times and reporting multiple annotations to NCBI RefSeq for the
same site instance), we introduce meta-site concept. A meta-site is the collection of
all evidences for the same genomic binding location, distributed across multiple
reported sites that are a few bp off each other. Each meta-site is represented by one
of its members (leader site).

Evidence from two motif-associated sites is combined into one meta-site if the
overlap between two sites is larger than 75% of the combined site length.  A
non-motif associated site is integrated into a meta-site if it fully overlaps with
any motif-associated site in the meta-site.

Meta-sites are not stored in the database and generated dynamically. For NCBI
annotation, overlapping site instances are merged into a meta-site and leader-site id
is used for dbxref link-out.

dbxref ID generation
--------------------
Each CollecTF annotation in NCBI RefSeq records has unique dbxref out-link to
CollecTF record, pointing a curated site instance. There is one-to-one mapping
between each annotated site-instance id and corresponding dbxref.
