from django.db import models
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
from django.core.cache import cache
import sys
import bioutils

# Create your models here.
class Curation(models.Model):
    """This table contains all the details about the curation, such as reported
    TF and species, followed experimental process, link to curator, publication,
    etc. It also keeps some meta-information about the curation, such as whether
    it requires revision, ready for NCBI submission, validation status, etc."""
    # choices
    REVISION_REASONS = (("genome_not_available", "No comparable genome in NCBI"),
                        ("in_progress", "Matching genome still in progress"),
                        ("TF_not_available", "No comparable TF protein sequence in NCBI"),
                        ("external_submission", "External submission"),
                        ("other", "other reason (specify in notes)"),)

    curation_id = models.AutoField(primary_key=True)

    # curation fields
    TF_species = models.CharField(max_length=500)   # species of reported TF
    site_species = models.CharField(max_length=500) # species of reported sites
    experimental_process = models.TextField(null=True, blank=True)

    # Is the TF shown to interact with other protein/ligand that influences
    # binding?
    forms_complex = models.BooleanField()
    complex_notes = models.TextField(null=True, blank=True)
    
    # curation meta information
    notes = models.TextField(blank=True)
    confidence = models.BooleanField()              # is curation confident?
    NCBI_submission_ready = models.BooleanField()   # is ready to submit to NCBI?
    requires_revision = models.CharField(max_length=20, choices=REVISION_REASONS, null=True, blank=True)
    validated_by = models.ForeignKey("Curator", null=True, blank=True, related_name="validated_by")

    # time stamps
    created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now=True)

    # relations
    curator = models.ForeignKey("Curator")
    publication = models.ForeignKey("Publication")
    TF = models.ForeignKey("TF", null=True)

    # <Wed Jan 29 2014> 1:N relationship between TF-instance and curation will
    # be changed into a N:N relationship. This will allow to define TFs composed
    # of several subunits with different accession numbers (e.g. heterodimer
    # such as IHF, composed of IHF alpha and beta.)
    #TF_instance = models.ForeignKey("TFInstance", related_name='x')
    TF_instances = models.ManyToManyField("TFInstance")

    # to be removed
    #TF_function = models.CharField(max_length=50)
    #TF_type = models.CharField(max_length=50)
    
    
    # <Wed Jan 29 2014> Experimental techniques are now going to be associated
    # directly (N:N) with curation_site_instances, allowing a curation to
    # contain curation_site_instances that have different associated techniques.
    #experimental_techniques = models.ManyToManyField("ExperimentalTechnique", related_name='curation')
    
    site_instances = models.ManyToManyField("SiteInstance", through="Curation_SiteInstance")

    # ChIP link (NULL if curation is not from ChIP paper
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
    """Curator table. Links to Django user model that contains basic stuff like
    name, email and password."""

    CURATOR_TYPE = (("internal", "internal"),
                    ("external", "external"))
    
    curator_id = models.AutoField(primary_key=True)
    user = models.OneToOneField(User) # extend Django's user model
    curator_type = models.CharField(max_length=20, choices=CURATOR_TYPE, default="external")

    def __unicode__(self):
        return u'%s' % self.user

class Publication(models.Model):
    """Publication table. It keeps publication information such as PMID,
    authors, title etc. It also has fields on whether it has promoter/expression
    data and reported TF/species that might be useful when the paper is
    curated"""
    
    PUBLICATION_TYPE = (("pubmed", "Pubmed article"),
                        ("nonpubmed", "Non-pubmed article"),
                        ("nonpublished", "Non-published data"))
    
    publication_id = models.AutoField(primary_key=True)
    publication_type = models.CharField(max_length=20, choices=PUBLICATION_TYPE)
    pmid = models.CharField(max_length=30, null=True, blank=True) # null if not PubMed article
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
    """Gene table containing position/strand on the genome, locus tag, accession
    number, etc."""
    
    STRAND = ((1, "Top strand"),
              (-1, "Bottom strand"))
    
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
    """Genome table"""
    GENOME_TYPE = (("chromosome", "chromosome"),
                   ("plasmid", "plasmid"))
    
    genome_id = models.AutoField(primary_key=True)
    genome_accession = models.CharField(max_length=20, unique=True)
    genome_sequence = models.OneToOneField("GenomeSequence", null=False, blank=False)
    GC_content = models.FloatField()
    gi = models.CharField(max_length=50, null=False)
    chromosome = models.CharField(max_length=10, null=False)
    organism = models.CharField(max_length=500, null=False)
    taxonomy = models.ForeignKey('Taxonomy')
    
    def __unicode__(self):
        return self.genome_accession + ' ' + self.organism

    def get_sequence(self):
        """Get genome sequence from the cache."""
        key = "genome-sequence-%s" % self.genome_accession
        if not cache.has_key(key):
            print "Key %s not in cache, retrieving from the database." % key
            value = self.genome_sequence.sequence
            value = str(value) # no need for unicode, less memory usage
            cache.set(key, value)

            cache.set('testing', '3')
            print cache.get('testing')
        ret = cache.get(key)
        assert ret
        return ret

    def get_genes(self):
        """Get the list of genes from the cache."""
        key = "genome-genes-%s" % self.genome_accession
        if not cache.has_key(key):
            print "Key %s not in cache, retrieving from the database." % key
            value = Gene.objects.filter(genome=self).order_by('start')
            cache.set(key,value)
        ret = cache.get(key)
        assert ret
        return ret
        

class GenomeSequence(models.Model):
    """Genome sequence table. It contains only sequence. Initially it was a part
    of genome table, but moved to a separate table for fast access to the genome
    table when the sequence is not necessary"""
    sequence = models.TextField(editable=False)
    
    def __unicode__(self):
        return '%s' % self.genome

class Taxonomy(models.Model):
    """Taxonomy table stores the phylogeny in the database."""
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
    
    def get_order(self):
        order = None
        x = self
        while (x.rank != 'order' and x.parent):
            x = x.parent
        if x.rank=='order':
            return x

    class Meta:
        verbose_name_plural = 'taxonomies'

class TF(models.Model):
    """TF table containing information about TF and link to its family"""
    TF_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    family = models.ForeignKey("TFFamily")
    description = models.TextField()

    def __unicode__(self):
        return u'%s [family: %s]' % (self.name, self.family.name)
    
    class Meta:
        verbose_name_plural = "TFs"

class TFFamily(models.Model):
    """TF family contains name and description for the family"""
    TF_family_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __unicode__(self):
        return u'%s' % self.name
    
    class Meta:
        verbose_name = "TF family"
        verbose_name_plural = "TF families"

class TFInstance(models.Model):
    """Protein table. Contains accession number, which TF it is and the
    description"""
    protein_accession = models.CharField(max_length=20, primary_key=True)
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __unicode__(self):
        return u'%s -- %s' % (self.protein_accession, self.description)

    class Meta:
        verbose_name = "TF instance"
        
class SiteInstance(models.Model):
    """The binding site model. Contains position/strand, which genome it is in
    etc."""
    site_id = models.AutoField(primary_key=True)
    _seq = models.TextField(max_length=100000) # redundant info, kept for sanity check
    genome = models.ForeignKey("Genome")
    # genome start position (0 index)
    start = models.IntegerField()
    # genome end position (end position! Not the first position following site
    # sequence, 0 index too.)
    end = models.IntegerField()
    # genome strand (1 or -1)
    strand = models.IntegerField(choices=Gene.STRAND) 

    def __unicode__(self):
        return u'%s [%s]' % (self.site_id, self._seq)

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

    def get_genome_sequence(self):
        return self.genome.get_sequence()
    

    @property
    def seq(self):
        genome = self.get_genome_sequence()
        sequence = genome[self.start:self.end+1]
        if self.strand == -1:
            # reverse complement
            sequence = bioutils.reverse_complement(sequence)
        assert sequence == self._seq
        return sequence

    @property
    def seq_lower(self):
        return str(self.seq).lower()
 
class Curation_SiteInstance(models.Model):
    """Through model between Curation and SiteInstance models. A specific
    binding site is unique (it is just start/end and strand afterall), but it
    might be determined in multiple papers. This table links curation and site
    instance tables and provides some additional information"""

    # TODO - explain site types.
    SITE_TYPE = (('motif_associated', "motif associated"),
                 ('non_motif_associated', "non-motif associated"),
                 ('var_motif_associated', "variable spaced motif associated"))

    TF_FUNCTION = (("ACT", "activator"),
                   ("REP", "repressor"),
                   ("DUAL", "dual"),
                   ("N/A", "not specified"))

    TF_TYPE = (("MONOMER", "monomer"),
               ("DIMER", "dimer"),
               ("TETRAMER", "tetramer"),
               ("OTHER", "other"),
               ("N/A", "not specified"))

    curation = models.ForeignKey("Curation", null=False)
    site_instance = models.ForeignKey("SiteInstance", null=False)
    site_type = models.CharField(max_length=50, choices=SITE_TYPE)
    annotated_seq = models.TextField(max_length=100000)

    # TF type
    TF_type = models.CharField(max_length=50, choices=TF_TYPE)
    # TF function
    TF_function = models.CharField(max_length=50, choices=TF_FUNCTION)

    # regulation relation
    regulates = models.ManyToManyField("Gene", through="Regulation")
    # experimental techniques used to determine the site instance
    experimental_techniques = models.ManyToManyField("ExperimentalTechnique")
    
    # If the paper reports sites that are identified by ChIP-based methods, each
    # site instance may have a quantitative value (e.g. peak intensity).
    quantitative_value = models.FloatField(null=True, blank=True)

    # If it is a high-throughput sequence, put a tag
    is_high_throughput = models.BooleanField(default=False)

    # NCBI submission related
    is_obsolete = models.BooleanField(null=False, default=False, blank=True)
    why_obsolete = models.TextField(null=True, blank=True) # explains why this site became obsolete.

    def __unicode__(self):
        return u'[%d] curation:%d' % (self.pk, self.curation.pk)

    @property
    def TF_function_verbose(self):
        return dict(Curation_SiteInstance.TF_FUNCTION)[self.TF_function]

    @property
    def TF_type_verbose(self):
        return dict(Curation_SiteInstance.TF_TYPE)[self.TF_type]

class Regulation(models.Model):
    """This table stores the TF regulation of genes. The regulation table links
    two tables: curation_site_instance and gene, to capture the information of
    which TF regulates which gene by binding which site on the genome. In
    addition, the evidence type is stored. If there is experimental evidence of
    regulation in the paper, saying that TF x up/down-regulates gene y, the
    evidence type is "exp-verified". Otherwise, all genes in the operon of which
    the site is upstream, are labeled as "inferred"."""
    
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
    """For some curations, the binding site that is reported in the paper can not be
    matched to any sequence in the reference genome for some reason (e.g. the
    reported genome is not available in RefSeq). Therefore, it is stored as
    not-annotated site instance which basically holds the reported sequence and
    a link to the curation table."""
    sequence = models.TextField(max_length=100000)
    curation = models.ForeignKey("Curation")
    
    def __unicode__(self):
        return u'%s [%s]' % (self.id, self.sequence)  
    
class ExperimentalTechnique(models.Model):
    """The table of experimental techniques. Each binding site is identified
    experimentally using a set of techniques. Curation of each site should have
    the information that how that particular site is experimentally
    determined. This table contains all the experimental techniques that can be
    used for TFBS identification, along with their description, type (used to
    show binding/expression or just in-silico)."""
    
    FUNCTION_CATEGORIES = (("binding", "Detection of binding"),
                           ("expression", "Assessment of expression"),
                           ("insilico", "In-silico prediction"))
    
    technique_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    description = models.TextField()
    preset_function = models.CharField(max_length=50, choices=FUNCTION_CATEGORIES, null=True)
    categories = models.ManyToManyField("ExperimentalTechniqueCategory")

    def __unicode__(self):
        return u'%s' % self.name

class ExperimentalTechniqueCategory(models.Model):
    """Each experimental technique can belong to a set of categories. Categories
    can be considered as labels that define the technique. It is different from
    preset_function field in the experimental_technique table which can have
    only three different values (binding, expression and in-silico), because an
    experimental technique may have several category. For instance, the
    technique ChIP-PCR has three categories: (1) PCR-based techniques,
    (2) Immunoprecipitation methods and (3) in vivo techniques. """
    category_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    description = models.TextField()

    def __unicode__(self):
        return u'[%d] %s' % (self.category_id, self.name)
    
    class Meta:
        verbose_name_plural = "experimental technique categories"

class ChipInfo(models.Model):
    """Some papers report binding sites using ChIP techniques. For these
    curations, this table contains information on ChIP experiment details. The
    curation table has a link (chip_info) to this table and it is NULL if the
    paper does not have ChIP technique."""
    chip_info_id = models.AutoField(primary_key=True)
    assay_conditions = models.TextField()
    method_notes = models.TextField()

    def __unicode__(self):
        return u'[%d] %s' % (self.chip_info_id, self.assay_conditions[:20])

class ExternalDatabase(models.Model):
    """Sometimes, authors choose to upload additional data (e.g. DNA-array data
    or gene expression data) to a database. To capture that information, on the
    submission process the curator is asked to provide the external database
    name and the accession number of the specific data reported. This table
    stores the information to access reported additional data."""

    ext_database_id = models.AutoField(primary_key=True)
    ext_database_name = models.CharField(max_length=50, null=False, unique=True)
    ext_database_descripton = models.CharField(max_length=500)
    ext_database_url_format = models.CharField(max_length=500)
    
    def __unicode__(self):
        return u'%s' % self.ext_database_name

class Curation_ExternalDatabase(models.Model):
    """This is just an intermediate table between curation table and the
    external database table. It links curation table to the external DB table
    and stores the accession number to the reported data."""
    curation = models.ForeignKey(Curation)
    external_database = models.ForeignKey(ExternalDatabase)
    accession_number = models.CharField(max_length=500, null=False)
    
    def __unicode__(self):
        return u'curation: %d - xref: %s [%s]' % (self.curation.curation_id,
                                                  self.external_database.ext_database_name,
                                                  self.accession_number)
class NCBISubmission(models.Model):
    """The curated data in CollecTF is integrated into NCBI RefSeq database
       through periodic genome-specific submissions. This internal table keeps
       track of all the data that has been submitted to the NCBI RefSeq."""

    submission_time = models.DateTimeField(auto_now_add=True)
    genome_submitted_to = models.CharField(max_length=50)
    curation_site_instance = models.ForeignKey('Curation_SiteInstance')
    
    class Meta:
        verbose_name = 'NCBI Submission'
