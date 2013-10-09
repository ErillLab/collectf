from django.db import models
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
from baseapp import bioutils

# Create your models here.
class Curation(models.Model):
    """
    Curation table.
    """
    # choices
    REVISION_REASONS = (("genome_not_available", "No comparable genome in NCBI"),
                        ("in_progress", "Matching genome still in progress"),
                        ("TF_not_available", "No comparable TF protein sequence in NCBI"),
                        ("external_submission", "External submission"),
                        ("other", "other reason (specify in notes)"),)    #other reasons
    
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
    TF_species = models.CharField(max_length=500)   # species of reported TF
    site_species = models.CharField(max_length=500) # species of reported sites
    confidence = models.BooleanField()              # is curation confident?
    NCBI_submission_ready = models.BooleanField()   # is ready to submit to NCBI?
    requires_revision = models.CharField(max_length=20, choices=REVISION_REASONS, null=True, blank=True)
    experimental_process = models.TextField(null=True, blank=True)
    forms_complex = models.BooleanField()           # does TF forms complex
    complex_notes = models.TextField(null=True, blank=True) # if forms complex,
    notes = models.TextField(blank=True)
    created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now=True)

    # Same TF can be both activator and repressor for different
    # situations. Similarly, same TF protein can bind DNA both as monomer and
    # dimer, etc. Therefore, instead of storing TF function and type in
    # TFInstance model, keep them here
    TF_function = models.CharField(max_length=50, choices=TF_FUNCTION)
    TF_type = models.CharField(max_length=50, choices=TF_TYPE)

    # relations
    curator = models.ForeignKey("Curator")
    publication = models.ForeignKey("Publication")
    TF = models.ForeignKey("TF", null=True)
    TF_instance = models.ForeignKey("TFInstance")
    experimental_techniques = models.ManyToManyField("ExperimentalTechnique", related_name='curation')
    site_instances = models.ManyToManyField("SiteInstance", through="Curation_SiteInstance")

    # ChIP link (NULL if site instance is not curated as ChIP data)
    chip_info = models.ForeignKey("ChipInfo", null=True, blank=True)
    quantitative_data_format = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return u'%s - %s - %s, %s, %s' % (self.curation_id,
                                          self.TF.name,
                                          self.publication.title,
                                          self.publication.authors,
                                          self.publication.publication_date)

    def _get_display(self, key, list):
        d = dict(list)
        if key in d: return d[key]
        return None
    
    def TF_function_verbose(self):
        return self._get_display(self.TF_function, self.TF_FUNCTION)

    def TF_type_verbose(self):
        return self._get_display(self.TF_type, self.TF_TYPE)
        
    

class Curator(models.Model):
    """
    Curator table. Each curator entry is linked to Django user object
    """
    curator_id = models.AutoField(primary_key=True)
    user = models.OneToOneField(User) # extend Django's user model
    is_staff = models.BooleanField(default=False, blank=True)

    def __unicode__(self):
        return u'%s' % self.user

class Publication(models.Model):
    """
    """
    PUBLICATION_TYPE = (("pubmed", "Pubmed article"),
                        ("nonpubmed", "Non-pubmed article"),
                        ("nonpublished", "Non-published data"))
    publication_id = models.AutoField(primary_key=True)
    publication_type = models.CharField(max_length=20, choices=PUBLICATION_TYPE)
    pmid = models.CharField(max_length=30, null=True, blank=True)
    # pmid can be null if not pubmed article
    authors = models.CharField(max_length=1000)
    title = models.CharField(max_length=1000)
    journal = models.CharField(max_length=1000)
    publication_date = models.CharField(max_length=50)
    volume = models.CharField(max_length=50)
    issue = models.CharField(max_length=50)
    pages = models.CharField(max_length=50)
    url = models.CharField(max_length=1000, null=True, blank=True)
    pdf = models.FileField(upload_to="papers/", null=True, blank=True)
    contains_promoter_data = models.BooleanField()
    contains_expression_data = models.BooleanField()
    submission_notes = models.TextField(null=True, blank=True)
    curation_complete = models.BooleanField(default=False) #paper curated?
    assigned_to = models.ForeignKey("Curator", null=True, blank=True)
    reported_TF = models.CharField(max_length=100, null=True, blank=True)
    reported_species = models.CharField(max_length=100, null=True, blank=True)
    
    def __unicode__(self):
        return u'[%s] PMID: %s, TF: %s, species: %s,assigned to: %s' % \
               (self.publication_id, self.pmid, self.reported_TF,
                self.reported_species, self.assigned_to)

class Gene(models.Model):
    """
    """
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
    def __unicode__(self):
        return '%s (%s-%s)' % (self.gene_id, self.name, self.genome.genome_accession)

class Genome(models.Model):
    """
    """
    GENOME_TYPE = (("chromosome", "chromosome"), ("plasmid", "plasmid"))
    genome_id = models.AutoField(primary_key=True)

    """from S.Lott@stackoverflow: The formal primary key should always be a surrogate
    key. Never anything else. [Strong words. Been database designer since the
    1980's. Important lessoned learned is this: everything is changeable, even when
    the user's swear on their mother's graves that the value cannot be changed is is
    truly a natural key that can be taken as primary. It isn't primary. Only
    surrogates can be primary.]"""
    genome_accession = models.CharField(max_length=20, unique=True)
    sequence = models.TextField(editable=False)
    GC_content = models.FloatField()
    gi = models.CharField(max_length=50, null=False)
    chromosome = models.CharField(max_length=10, null=False)
    organism = models.CharField(max_length=500, null=False)
    taxonomy = models.ForeignKey('Taxonomy')
    
    def __unicode__(self):
        return self.genome_accession + ' ' + self.organism

class Taxonomy(models.Model):
    """
    """
    taxonomy_id = models.CharField(max_length=20, unique=True)
    rank = models.CharField(max_length=20, choices=(('phylum', 'phylum'),
                                                    ('class', 'class'),
                                                    ('order', 'order'),
                                                    ('family', 'family'),
                                                    ('genus', 'genus'),
                                                    ('species', 'species')), null=False)
    name = models.CharField(max_length=100)
    parent = models.ForeignKey('self', null=True)
    def __unicode__(self):
        return '[%s] %s (%s)' % (str(self.taxonomy_id), self.name, self.rank)
    class Meta:
        verbose_name_plural = 'taxonomies'

    def get_order(self):
        order = None
        x = self
        while (x.rank != 'order' and x.parent):
            x = x.parent
        if x.rank=='order':
            return x
            

class TF(models.Model):
    """
    """
    TF_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    family = models.ForeignKey("TFFamily")
    description = models.TextField()

    def __unicode__(self):
        return u'%s (family: %s)' % (self.name, self.family.name)
    
    class Meta:
        verbose_name_plural = "TFs"

class TFFamily(models.Model):
    
    TF_family_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __unicode__(self):
        return u'%s' % self.name
    
    class Meta:
        verbose_name = "TF family"
        verbose_name_plural = "TF families"

class TFInstance(models.Model):
    protein_accession = models.CharField(max_length=20, primary_key=True)
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __unicode__(self):
        return u'%s -- %s' % (self.protein_accession, self.description)

    class Meta:
        verbose_name = "TF instance"
        
class SiteInstance(models.Model):
    site_id = models.AutoField(primary_key=True)
    seq = models.TextField(max_length=100000)
    genome = models.ForeignKey("Genome")
    start = models.IntegerField() # genome start position (0 index)
    end = models.IntegerField()   # genome end position (end position! not the first position after site sequence, 0 index too.)
    strand = models.IntegerField(choices=Gene.STRAND) # genome strand (1 or -1)

    def __unicode__(self):
        return u'%s [%s]' % (self.site_id, self.seq)

    def to_fasta(self):
        desc = "%s %s(%d, %d)" % (self.genome.genome_accession,
                                  '+' if self.strand==1 else '-',
                                  self.start,
                                  self.end)
        seq = self.seq
        return ">%s\n%s\n" % (desc, seq)

    def to_csv(self):
        fields = [self.genome.genome_accession,
                  '+' if self.strand == 1 else '-',
                  '%d' % self.start,
                  '%d' % self.end,
                  self.seq,]
        return '\t'.join(fields)

class Curation_SiteInstance(models.Model):
    # through model between Curation and SiteInstance models
    curation = models.ForeignKey("Curation", null=False)
    site_instance = models.ForeignKey("SiteInstance", null=False)
    is_motif_associated = models.BooleanField(null=False) # is the site instance actual motif,
                                                 # or a (longer) sequence that
                                                 # contains the site somewhere
    annotated_seq = models.TextField(max_length=100000)
    # regulation
    regulates = models.ManyToManyField("Gene", through="Regulation")
    # quantitative value
    quantitative_value = models.FloatField(null=True, blank=True)

    # NCBI submission related
    is_obsolete = models.BooleanField(null=False, default=False, blank=True)
    why_obsolete = models.TextField(null=True, blank=True) # explain why this site became obsolete.


    def __unicode__(self):
        return u'[%d]' % self.pk
        #return u"reported: %s, matched: %s" % (self.site_instance, self.annotated_seq)

class Regulation(models.Model):
    EVIDENCE_TYPE = (("exp_verified", "experimentally verified"),
                     ("inferred", "inferred"))
    curation_site_instance = models.ForeignKey("Curation_SiteInstance")
    gene = models.ForeignKey("Gene")
    evidence_type = models.CharField(max_length=20, choices=EVIDENCE_TYPE)

    def __unicode__(self):
        return 'curation_id: %s gene: %s, site_id: %s, type: %s' % \
               (self.curation_site_instance.curation.curation_id,
                self.gene.name, self.curation_site_instance.site_instance.site_id,
                self.evidence_type)

class NotAnnotatedSiteInstance(models.Model):
    """If no matching sequence found in genome, use this class"""
    sequence = models.TextField(max_length=100000)
    curation = models.ForeignKey("Curation")
    def __unicode__(self):
        return u'%s [%s]' % (self.id, self.sequence)  
    
class ExperimentalTechnique(models.Model):
    technique_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    description = models.TextField()
    FUNCTION_CATEGORIES = (("binding", "Detection of binding"),
                           ("expression", "Assessment of expression"),
                           ("insilico", "In-silico prediction"))
    preset_function = models.CharField(max_length=50, choices=FUNCTION_CATEGORIES, null=True)
    categories = models.ManyToManyField("ExperimentalTechniqueCategory")

    def __unicode__(self):
        return u'%s' % self.name

class ExperimentalTechniqueCategory(models.Model):
    category_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    description = models.TextField()

    def __unicode__(self):
        return u'[%d] %s' % (self.category_id, self.name)
    
    class Meta:
        verbose_name_plural = "experimental technique categories"

class ChipInfo(models.Model):
    chip_info_id = models.AutoField(primary_key=True)
    assay_conditions = models.TextField()
    method_notes = models.TextField()

    def __unicode__(self):
        return u'[%d] PMID: %s' % (self.chip_info_id, self.curation_set.all()[0].publication.pmid)

class ExternalDatabase(models.Model):
    ext_database_id = models.AutoField(primary_key=True)
    ext_database_name = models.CharField(max_length=50, null=False, unique=True)
    ext_database_descripton = models.CharField(max_length=500)
    ext_database_url_format = models.CharField(max_length=500)
    def __unicode__(self):
        return u'%s' % self.ext_database_name

class Curation_ExternalDatabase(models.Model):
    curation = models.ForeignKey(Curation)
    external_database = models.ForeignKey(ExternalDatabase)
    accession_number = models.CharField(max_length=500, null=False)
    def __unicode__(self):
        return u'curation: %d - xref: %s [%s]' % (self.curation.curation_id,
                                                  self.external_database.ext_database_name,
                                                  self.accession_number)
class NCBISubmission(models.Model):
    submission_time = models.DateTimeField(auto_now_add=True)
    genome_submitted_to = models.CharField(max_length=50)
    curation_site_instance = models.ForeignKey('Curation_SiteInstance')
    class Meta:
        verbose_name = 'NCBI Submission'
