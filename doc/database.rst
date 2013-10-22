Database Structure
==================

Graphical representation of the database is available at :download:`collectf_db.pdf`.

The database consists of several tables and relationships between them. It is
implemented in MySQL with Django interface. Each table is represented as a Python
class (model) and access to the database is through instances of those classes.

Below, each table is described in detail.

Curation
--------

Curator
-------

Publication
-----------

Gene
----

Genome
------

GenomeSequence
--------------

Taxonomy
--------

TF
--

TFFamily
--------

TFInstance
----------

SiteInstance
------------

Curation_SiteInstance
---------------------

Regulation
----------

NotAnnotatedSiteInstance
------------------------

ExperimentalTechnique
---------------------

ChipInfo
--------

ExternalDatabase
----------------

Curation_ExternalDatabase
-------------------------

NCBISubmission
--------------

