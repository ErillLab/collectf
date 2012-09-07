from django.db import models
from django.contrib.auth.models import User

# Create your models here.
class Curation(models.Model):
    # choices
    REVISION_REASONS = (("not_available", "genome not available"),
                        ("not_sequenced", "sequence not available"),
                        ("in_progress", "in progress"), #genome is being sequenced
                        ("other", "other"),)            #other reasons
    
    TF_FUNCTION = (("ACT", "activator"),
                   ("REP", "repressor"),
                   ("N/A", "not specified"))
    
    TF_TYPE = (("MONOMER", "monomer"),
               ("DIMER", "dimer"),
               ("TETRAMER", "tetramer"),
               ("OTHER", "other"),
               ("N/A", "not specified"))
    # fields
    curation_id = models.AutoField(primary_key=True)
    TF_species = models.CharField(max_length=100)   # species of reported TF
    site_species = models.CharField(max_length=100) # species of reported sites
    confidence = models.BooleanField()              # is curation confident?
    NCBI_submission_ready = models.BooleanField()   # is ready to submit to NCBI?
    requires_revision = models.CharField(max_length=20, choices=REVISION_REASONS,
                                         null=True)
    experimental_process = models.TextField(null=True, blank=True)
    forms_complex = models.BooleanField()           # does TF forms complex
    complex_notes = models.TextField(null=True, blank=True) # if forms complex,
    notes = models.TextField()
    last_modified = models.DateTimeField()

    # Same TF can be both activator and repressor for different
    # situations. Similarly, same TF protein can bind DNA both as monomer and
    # dimer, etc. Therefore, instead of storing TF function and type in
    # TFInstance model, keep them here
    TF_function = models.CharField(max_length=5, choices=TF_FUNCTION)
    TF_type = models.CharField(max_length=5, choices=TF_TYPE)

    # relations
    curator = models.ForeignKey("Curator")
    publication = models.ForeignKey("Publication")
    TF = models.ForeignKey("TF", null=True)
    TF_instance = models.ForeignKey("TFInstance")
    experimental_techniques = models.ManyToManyField("ExperimentalTechnique")
    site_instances = models.ManyToManyField("SiteInstance")
    
class Curator(models.Model):
    curator_id = models.AutoField(primary_key=True)
    user = models.OneToOneField(User) # extend Django's user model
    assigned_papers = models.ManyToManyField("Publication")

class Publication(models.Model):
    PUBLICATION_TYPE = (("pubmed", "Pubmed article"),
                        ("nonpubmed", "Non-pubmed article"),
                        ("nonpublished", "Non-published data"))
    publication_id = models.AutoField(primary_key=True)
    publication_type = models.CharField(max_length=20, choices=PUBLICATION_TYPE)
    pmid = models.CharField(max_length=30, null=True, blank=True)
    # pmid can be null if not pubmed article
    citation = models.TextField() # contains authors, title and journal info
    url = models.URLField(null=True, blank=True)
    pdf = models.FileField(upload_to="papers/", null=True, blank=True)
    contains_promoter_data = models.BooleanField()
    contains_expression_data = models.BooleanField()
    submission_notes = models.TextField(null=True, blank=True)
    curation_complete = models.BooleanField(default=False) #paper curated?

class Gene(models.Model):
    # choices
    STRAND = ((1, "Top strand"), (-1, "Bottom strand"))
    # fields
    gene_id = models.AutoField(primary_key=True)
    gene_accession = models.CharField(max_length=30, null=True, blank=True)
    genome = models.ForeignKey("Genome")
    name = models.CharField(max_length=50)
    description = models.CharField(max_length=1000)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField(choices=STRAND)
    locus_tag = models.CharField(max_length=20)
    # COG id?
    # KEGG id?
    # homology

class Genome(models.Model):
    GENOME_TYPE = (("chromosome", "chromosome"), ("plasmid", "plasmid"))
    genome_accession = models.CharField(max_length=20, primary_key=True)
    genome_type = models.CharField(max_length=20, choices=GENOME_TYPE)
    sequence = models.TextField()
    GC_content = models.FloatField()
    strain = models.ForeignKey("Strain")

class Strain(models.Model):
    taxonomy_id = models.CharField(max_length=200, primary_key=True)
    name = models.CharField(max_length=100)

class TF(models.Model):
    TF_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    family = models.ForeignKey("TFFamily")
    description = models.TextField()
    class Meta:
        verbose_name_plural = "TFs"

class TFFamily(models.Model):
    TF_family_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    description = models.TextField()
    class Meta:
        verbose_name = "TF family"
        verbose_name_plural = "TF families"

class TFInstance(models.Model):
    protein_accession = models.CharField(max_length=20, primary_key=True)
    name = models.CharField(max_length=50)
    description = models.TextField()

    class Meta:
        verbose_name = "TF instance"

class SiteInstance(models.Model):
    site_id = models.AutoField(primary_key=True)
    annotated_seq = models.TextField()
    original_seq = models.TextField()
    genome = models.ForeignKey("Genome")
    start = models.IntegerField() # genome start position
    end = models.IntegerField()   # genome end position
    strand = models.IntegerField(choices=Gene.STRAND) # genome strand (1 or -1)
    regulates = models.ManyToManyField("Gene", through="Regulation")

class Regulation(models.Model):
    EVIDENCE_TYPE = (("exp", "experimentally verified"), ("inferred", "inferred"))
    site_instance = models.ForeignKey("SiteInstance")
    gene = models.ForeignKey("Gene")
    evidence_type = models.CharField(max_length=20, choices=EVIDENCE_TYPE)

class NotAnnotatedSiteInstance(models.Model):
    """If no matching sequence found in genome, use this class"""
    sequence = models.TextField()
    curation = models.ForeignKey("Curation")
    
class ExperimentalTechnique(models.Model):
    technique_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    description = models.TextField()
    

