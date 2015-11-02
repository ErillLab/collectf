"""CollecTF Django model definitions."""

import Queue

from django.db import models
from django.contrib.auth.models import User
from django.core.cache import cache

class Curation(models.Model):
    """Curation model.

    Contains all the details about the curation, such as reported TF and
    species, followed experimental process, link to curator, publication,
    etc. Also keeps some meta-information about the curation, such as whether it
    requires revision, ready for NCBI submission, validation status, etc.
    """
    # Curation identifier
    curation_id = models.AutoField(primary_key=True)

    # Species of the report transcription factor
    TF_species = models.CharField(max_length=500)

    # Species of the report species
    site_species = models.CharField(max_length=500)

    # Description of the experimental process
    experimental_process = models.TextField(null=True, blank=True)

    # If true, TF is shown to interact with other protein/ligand.
    forms_complex = models.BooleanField()

    # Description of the interaction with other protein/ligand.
    complex_notes = models.TextField(null=True, blank=True)

    # Curation notes
    notes = models.TextField(blank=True)

    # If true, the curator is confident about the curation.
    confidence = models.BooleanField()

    # If true, the curation is ready for NCBI submission.
    NCBI_submission_ready = models.BooleanField()

    # Reason for the revision of the curation.
    REVISION_REASONS = (
        ('genome_not_available', "No comparable genome in NCBI"),
        ('in_progress', "Matching genome still in progress"),
        ('TF_not_available', "No comparable TF protein sequence in NCBI"),
        ('other', "Other reason (specify in notes)"),)
    requires_revision = models.CharField(
        max_length=20, choices=REVISION_REASONS, null=True, blank=True)

    # The curator that validated the curation.
    validated_by = models.ForeignKey('Curator', null=True, blank=True,
                                     related_name='validated_by')

    # The time that the curation is created.
    created = models.DateTimeField(auto_now_add=True)

    # The time that the curation is last modified.
    last_modified = models.DateTimeField(auto_now=True)

    # The curator that submitted the curation.
    curator = models.ForeignKey('Curator')

    # The publication that the curation belongs to.
    publication = models.ForeignKey('Publication')

    # The transcription factor instances that the curation has.
    TF_instances = models.ManyToManyField('TFInstance')

    # The binding sites
    site_instances = models.ManyToManyField('SiteInstance',
                                            through='Curation_SiteInstance')

    # ChIP experiment data if the curation belongs to a ChIP paper.
    chip_info = models.ForeignKey('ChipInfo', null=True, blank=True)

    # The format of the quantitative data the curation has.
    quantitative_data_format = models.TextField(null=True, blank=True)

    def __unicode__(self):
        """Returns a unicode representation of the curation object."""
        return u'%s - %s - %s, %s, %s' % (self.curation_id,
                                          self.TF.name,
                                          self.publication.title,
                                          self.publication.authors,
                                          self.publication.publication_date)

    @property
    def TF(self):
        """Returns the TF associated with the curation."""
        return self.TF_instances.all()[0].TF

    def TF_function_verbose(self):
        return dict(self.TF_FUNCTION)[self.TF_function]

    def TF_type_verbose(self):
        return dict(self.TF_TYPE)[self.TF_type]

class Curator(models.Model):
    """Curator table.

    Links to Django user model that contains name, email and password.
    """
    # The curator identifier
    curator_id = models.AutoField(primary_key=True)

    # The Django user object that the curator is linked to.
    user = models.OneToOneField(User)

    # The type of the curator
    CURATOR_TYPE = (('internal', 'internal'), ('external', 'external'))
    curator_type = models.CharField(max_length=20, choices=CURATOR_TYPE,
                                    default='external')

    def __unicode__(self):
        """Returns the unicode representation of the curator."""
        return u'%s' % self.user

class Publication(models.Model):
    """Publication table."""

    PUBLICATION_TYPE = (('pubmed', "Pubmed article"),
                        ('nonpubmed', "Non-pubmed article"),
                        ('nonpublished', "Non-published data"))

    # Publication identifier.
    publication_id = models.AutoField(primary_key=True)

    # The type of the publication.
    publication_type = models.CharField(max_length=20, choices=PUBLICATION_TYPE)

    # The PubMed identifier of the publication.
    pmid = models.CharField(max_length=30, null=True, blank=True)

    # The authors of the publication.
    authors = models.CharField(max_length=1000)

    # The title of the publication.
    title = models.CharField(max_length=1000)

    # The journal that the paper is published.
    journal = models.CharField(max_length=1000)

    # The publication date.
    publication_date = models.CharField(max_length=50)

    # The volume, issue and pages of the publication.
    volume = models.CharField(max_length=50)
    issue = models.CharField(max_length=50)
    pages = models.CharField(max_length=50)

    # The URL for the publication
    url = models.CharField(max_length=1000, null=True, blank=True)

    # TODO(sefa): remove unused pdf field.
    pdf = models.FileField(upload_to='papers/', null=True, blank=True)

    # If true, the paper contains promoter data.
    contains_promoter_data = models.BooleanField()

    # If true, the paper contains expression data.
    contains_expression_data = models.BooleanField()

    # The notes about the submission of the paper to CollecTF.
    submission_notes = models.TextField(null=True, blank=True)

    # If true, the curation(s) for this paper is complete.
    curation_complete = models.BooleanField(default=False)

    # The curator that the paper is assigned to.
    assigned_to = models.ForeignKey('Curator', null=True, blank=True)

    # The transcription factor that is reported in the paper.
    reported_TF = models.CharField(max_length=100, null=True, blank=True)

    # The species that is reported in the paper.
    reported_species = models.CharField(max_length=100, null=True, blank=True)

    def __unicode__(self):
        """Returns the unicode representation of the publication."""
        return (u'[%s] PMID: %s, TF: %s, species: %s,assigned to: %s' %
                (self.publication_id, self.pmid, self.reported_TF,
                 self.reported_species, self.assigned_to))

class Gene(models.Model):
    """Gene table.

    Contains position/strand on the genome, locus tag, accession number.
    """

    # The gene identifier.
    gene_id = models.AutoField(primary_key=True)

    # The NCBI accession number of the gene.
    gene_accession = models.CharField(max_length=30, null=True, blank=True)

    # The genome of the gene.
    genome = models.ForeignKey('Genome')

    # The name of the gene.
    name = models.CharField(max_length=50)

    # The description of the gene.
    description = models.CharField(max_length=1000)

    # The start position.
    start = models.IntegerField()

    # The end position.
    end = models.IntegerField()

    # The strand that the gene lies on.
    STRAND = ((1, "Top strand"), (-1, "Bottom strand"))
    strand = models.IntegerField(choices=STRAND)

    # The locus tag of the gene.
    locus_tag = models.CharField(max_length=20)

    def __unicode__(self):
        """Returns the unicode representation of the gene."""
        return '%s (%s-%s)' % (self.gene_id, self.name,
                               self.genome.genome_accession)

class Genome(models.Model):
    """Genome table."""
    # The genome identifier.
    genome_id = models.AutoField(primary_key=True)

    # The NCBI accession number of the genome.
    genome_accession = models.CharField(max_length=20, unique=True)

    # The genome sequence.
    genome_sequence = models.OneToOneField('GenomeSequence', null=False,
                                           blank=False)

    # The GC content of the genome sequence.
    GC_content = models.FloatField()

    # The GI number of the genome.
    gi = models.CharField(max_length=50, null=False)

    # The chromosome of the genome.
    chromosome = models.CharField(max_length=10, null=False)

    # The organism of the genome.
    organism = models.CharField(max_length=500, null=False)

    # The taxonomy of the species that this genome belongs to.
    taxonomy = models.ForeignKey('Taxonomy')

    def __unicode__(self):
        """Returns the unicode representation of the Genome."""
        return self.genome_accession + ' ' + self.organism

    def get_sequence(self):
        """Gets genome sequence from the cache."""
        key = 'genome_sequence_%s' % self.genome_accession
        if not cache.has_key(key):
            value = self.genome_sequence.sequence
            value = str(value) # No need for unicode, less memory usage.
            cache.set(key, value)
        ret = cache.get(key)
        return ret

    def get_genes(self):
        """Gets the list of genes from the cache."""
        key = 'genome_genes_%s' % self.genome_accession
        if not cache.has_key(key):
            value = Gene.objects.filter(genome=self).order_by('start')
            cache.set(key, value)
        ret = cache.get(key)
        return ret

class GenomeSequence(models.Model):
    """Genome sequence table.

    Contains only sequence. Was a part of genome table, but moved to a separate
    table for efficient access to the genome table when the sequence is not
    necessary.
    """

    # The genome sequence.
    sequence = models.TextField(editable=False)

    def __unicode__(self):
        """Unicode representation of the genome sequence."""
        return '%s' % self.genome

class Taxonomy(models.Model):
    """The phylogeny in the database."""

    # The taxonomy identifier.
    taxonomy_id = models.CharField(max_length=20, unique=True)

    # The rank of the taxonomy object.
    rank = models.CharField(max_length=20,
                            choices=(('phylum', 'phylum'),
                                     ('class', 'class'),
                                     ('order', 'order'),
                                     ('family', 'family'),
                                     ('genus', 'genus'),
                                     ('species', 'species')),
                            null=False)

    # The name of the taxonomy object.
    name = models.CharField(max_length=100)

    # The parent node.
    parent = models.ForeignKey('self', null=True)

    def __unicode__(self):
        """Returns the unicode representation of the taxonomy. object."""
        return '[%s] %s (%s)' % (str(self.taxonomy_id), self.name, self.rank)

    def get_all_species(self):
        """Returns all species of the subtree."""
        all_species = []
        node_queue = Queue.Queue()
        node_queue.put(self)
        while not node_queue.empty():
            node = node_queue.get()
            children = node.taxonomy_set.all()
            if children:
                for child in children:
                    node_queue.put(child)
            else:
                all_species.append(node)
        return all_species

    class Meta:
        verbose_name_plural = 'taxonomies'
        ordering = ['name']

class TF(models.Model):
    """Transcription factor and link to its family."""

    # The TF identifier.
    TF_id = models.AutoField(primary_key=True)

    # The name of the TF.
    name = models.CharField(max_length=50)

    # The TF family.
    family = models.ForeignKey('TFFamily')

    # The description of the TF.
    description = models.TextField()

    def __unicode__(self):
        """Returns the unicode representation of the TF."""
        return u'%s [family: %s]' % (self.name, self.family.name)

    class Meta:
        verbose_name_plural = "TFs"
        ordering = ['name']

class TFFamily(models.Model):
    """TF family table."""

    # The TF family identifier.
    TF_family_id = models.AutoField(primary_key=True)

    # The name of the TF family.
    name = models.CharField(max_length=50)

    # The description of the TF family.
    description = models.TextField()

    def __unicode__(self):
        """Returns the unicode representation of the TF family object."""
        return u'%s' % self.name

    class Meta:
        verbose_name = "TF family"
        verbose_name_plural = "TF families"
        ordering = ['name']

class TFInstance(models.Model):
    """TF instance (protein) table."""

    # The TF instance identifier.
    TF_instance_id = models.AutoField(primary_key=True)

    # The RefSeq accession number.
    refseq_accession = models.CharField(max_length=20)

    # The UniProt accession number.
    uniprot_accession = models.CharField(max_length=20, unique=True)

    # The description of the TF instance.
    description = models.TextField()

    # Transcription factor that the protein is.
    TF = models.ForeignKey('TF', null=False)

    # The notes -- mostly about RefSeq to UniProt migration.
    notes = models.TextField()

    def __unicode__(self):
        """Returns the unicode representation of the TF instance."""
        return u'%s (%s): %s' % (
            self.uniprot_accession, self.refseq_accession, self.description)

    class Meta:
        ordering = ['uniprot_accession']
        verbose_name = "TF instance"

class SiteInstance(models.Model):
    """The binding site table."""

    # The binding site identifier.
    site_id = models.AutoField(primary_key=True)

    # The binding site sequence. Not used, stored for sanity check.
    _seq = models.TextField(max_length=100000)

    # The genome that the binding site belongs to.
    genome = models.ForeignKey("Genome")

    # The genome start position. 0-indexed
    start = models.IntegerField()
    # The genome end position (i.e. position of last base of the binding
    # site). 0-indexed.
    end = models.IntegerField()

    # The genome strand. 1 or -1
    strand = models.IntegerField(choices=Gene.STRAND)

    def __unicode__(self):
        """Returns the unicode representation of binding site."""
        return u'%s [%s]' % (self.site_id, self._seq)

    def get_genome_sequence(self):
        """Returns the genome sequence that binding site belongs to."""
        return self.genome.get_sequence()

    @property
    def seq(self):
        """Returns the sequence of the binding site."""
        genome = self.get_genome_sequence()
        sequence = genome[self.start:self.end+1]
        if self.strand == -1:
            sequence = bioutils.reverse_complement(sequence)
        assert sequence == self._seq
        return sequence

    @property
    def seq_lower(self):
        """Returns the binding site in lower-case letters."""
        return str(self.seq).lower()

class Curation_SiteInstance(models.Model):
    """'Through' model between Curation and SiteInstance models.

    A binding site is unique (defined as the start and end positions and
    strand), but it might be reported by multiple papers. This table links
    curation and site instance tables and provides some additional information.
    """

    # The Curation object.
    curation = models.ForeignKey("Curation", null=False)

    # The SiteInstance object.
    site_instance = models.ForeignKey("SiteInstance", null=False)

    # The type of the binding site.
    SITE_TYPE = (('motif_associated', "motif associated"),
                 ('var_motif_associated', "variable motif associated"),
                 ('non_motif_associated', "non-motif associated"))
    site_type = models.CharField(max_length=50, choices=SITE_TYPE)

    # The binding site sequence as reported in the paper. It could be a
    # different sequence if a surrogate genome is used or if the reported
    # sequence is not found in the genome.
    annotated_seq = models.TextField(max_length=100000)

    # The type of the transcription factor.
    TF_TYPE = (('MONOMER', "monomer"),
               ('DIMER', "dimer"),
               ('TETRAMER', "tetramer"),
               ('OTHER', "other"),
               ('N/A', "not specified"))
    TF_type = models.CharField(max_length=50, choices=TF_TYPE)

    # The function of the transcription factor.
    TF_FUNCTION = (('ACT', "activator"),
                   ('REP', "repressor"),
                   ('DUAL', "dual"),
                   ('N/A', "not specified"))
    TF_function = models.CharField(max_length=50, choices=TF_FUNCTION)

    # The Regulation objects describing the genes regulated via this binding
    # site.
    regulates = models.ManyToManyField("Gene", through='Regulation')

    # The experimental techniques used to determine this site instance.
    experimental_techniques = models.ManyToManyField('ExperimentalTechnique')

    # The associated quantitative value that the binding site could have if the
    # paper reports sites that are identified through high-throughput methods
    # (e.g. peak density).
    quantitative_value = models.FloatField(null=True, blank=True)

    # True if the binding site is a high-throughput sequence.
    is_high_throughput = models.BooleanField(default=False)

    # True, if the binding site is obsolete. A site is marked as obsolete if it
    # has errors and can not be deleted since it has been submitted to NCBI.
    is_obsolete = models.BooleanField(null=False, default=False, blank=True)

    # The reason for marking the binding site as obsolete.
    why_obsolete = models.TextField(null=True, blank=True)

    # The motif identifier. A TF could recognize more than one motif per
    # species. This field stores which motif that the binding site belongs to.
    motif_id = models.IntegerField(default=-1, null=False, blank=False)

    def __unicode__(self):
        """Returns the unicode representation of the Curation_SiteInstance."""
        return u'[%d] curation:%d TF:%s species:%s' % (
            self.pk, self.curation.pk, self.curation.TF,
            self.site_instance.genome)

    @property
    def TF_function_verbose(self):
        return dict(Curation_SiteInstance.TF_FUNCTION)[self.TF_function]

    @property
    def TF_type_verbose(self):
        return dict(Curation_SiteInstance.TF_TYPE)[self.TF_type]

    @property
    def genome_accession(self):
        """Returns the genome accession of the site instance."""
        return self.site_instance.genome.genome_accession



class Regulation(models.Model):
    """Gene regulation table.

    Stores the TF regulation of genes. Links two tables, Curation_SiteInstance
    and Gene, to capture the information of which TF regulates which gene by
    binding which site on the genome.

    Also stores the evidence type. If there is experimental evidence of
    regulation in the paper, saying that TF x up/down-regulates gene y, the
    evidence type is 'exp-verified'. Otherwise, all genes in the operon of which
    the site is upstream, are labeled as 'inferred'.
    """
    # The Curation_SiteInstance object.
    curation_site_instance = models.ForeignKey('Curation_SiteInstance')

    # The regulated gene.
    gene = models.ForeignKey("Gene")

    # The type of evidence for the regulation.
    EVIDENCE_TYPE = (('exp_verified', "experimentally verified"),
                     ('inferred', "inferred"))
    evidence_type = models.CharField(max_length=20, choices=EVIDENCE_TYPE)

    def __unicode__(self):
        """Returns the unicode representation of the Regulation."""
        return 'curation_id: %s gene: %s, site_id: %s, type: %s' % (
            self.curation_site_instance.curation.curation_id,
            self.gene.name,
            self.curation_site_instance.site_instance.site_id,
            self.evidence_type)

class NotAnnotatedSiteInstance(models.Model):
    """The table for not annotated site instances.

    For some curations, the binding site that is reported in the paper could not
    be matched to any sequence in the reference genome for some reason (e.g. the
    reported genome is not available in RefSeq). Therefore, it is stored as
    NotAnnotatedSiteInstance which basically holds the reported sequence and a
    link to the curation table.
    """
    # The binding site sequence
    sequence = models.TextField(max_length=100000)

    # The curation of the binding site.
    curation = models.ForeignKey("Curation")

    # The type of the TF.
    TF_type = models.CharField(max_length=50,
                               choices=Curation_SiteInstance.TF_TYPE,
                               null=True)
    # The function of the TF.
    TF_function = models.CharField(max_length=50,
                                   choices=Curation_SiteInstance.TF_FUNCTION,
                                   null=True)

    # Experimental techniques used to determine the site instance.
    experimental_techniques = models.ManyToManyField('ExperimentalTechnique')

    # The associated quantitative value that the binding site could have if the
    # paper reports sites that are identified through high-throughput methods
    # (e.g. peak density).
    quantitative_value = models.FloatField(null=True, blank=True)

    # True if the binding site is a high-throughput sequence.
    is_high_throughput = models.BooleanField(default=False)

    def __unicode__(self):
        """Returns the unicode representation of the NotAnnotatedSiteInstance"""
        return u'%s [%s]' % (self.id, self.sequence)

class ExperimentalTechnique(models.Model):
    """The table of experimental techniques.

    Contains all the experimental techniques that can be used for TFBS
    identification, along with their description and type.
    """
    # The experimental technique identifier.
    technique_id = models.AutoField(primary_key=True)

    # The name of the technique.
    name = models.CharField(max_length=100)

    # The description of the technique.
    description = models.TextField()

    # The default function of the experimental technique.
    FUNCTION_CATEGORIES = (("binding", "Detection of binding"),
                           ("expression", "Assessment of expression"),
                           ("insilico", "In-silico prediction"))
    preset_function = models.CharField(max_length=50,
                                       choices=FUNCTION_CATEGORIES, null=True)

    # The categories of the technique.
    categories = models.ManyToManyField("ExperimentalTechniqueCategory")

    # Evidence ontology term for the technique.
    EO_term = models.CharField(max_length=50, null=True, blank=True)

    def __unicode__(self):
        """Returns the unicode representation of the experimental technique."""
        return u'%s' % self.name

    class Meta:
        ordering = ['name']

class ExperimentalTechniqueCategory(models.Model):
    """Experimental technique category table.

    Categories can be considered as labels that define the technique. It is
    different from preset_function field in the experimental_technique table
    which can have only three different values (binding, expression and
    in-silico), because an experimental technique may have several
    categories.

    For instance, the technique ChIP-PCR has three categories: (1) PCR-based
    techniques, (2) Immunoprecipitation methods and (3) in vivo techniques.
    """

    # The category identifier.
    category_id = models.AutoField(primary_key=True)

    # The name of the technique category.
    name = models.CharField(max_length=100)

    # The description of the technique category.
    description = models.TextField()

    def __unicode__(self):
        """Returns the unicode representation of the technique category."""
        return u'[%d] %s' % (self.category_id, self.name)

    class Meta:
        ordering = ['name']
        verbose_name_plural = "experimental technique categories"

class ChipInfo(models.Model):
    """Chromatin immunoprecipitation (ChIP) experiment information.

    Some papers report binding sites using ChIP techniques. For these curations,
    this table contains information on ChIP experiment details. The curation
    table has a link (chip_info) to this table and it is NULL if the paper does
    not have any ChIP technique.
    """

    # The identifier.
    chip_info_id = models.AutoField(primary_key=True)

    # The description of the assay conditions.
    assay_conditions = models.TextField()

    # The method notes.
    method_notes = models.TextField()

    def __unicode__(self):
        """Returns the unicode representation of the object."""
        return u'[%d] %s' % (self.chip_info_id, self.assay_conditions[:20])

class ExternalDatabase(models.Model):
    """The external database table.

    Sometimes, authors prefers uploading additional data (e.g. DNA-array data,
    gene expression data) to a database. To capture that information, during the
    submission process, the curator is asked to provide the external database
    name and the accession number of the specific data reported. This table
    stores the information to access reported additional data.
    """
    # External database identifier.
    ext_database_id = models.AutoField(primary_key=True)

    # The external database name.
    ext_database_name = models.CharField(max_length=50, null=False, unique=True)

    # The external database description.
    ext_database_descripton = models.CharField(max_length=500)

    # The external database URL format.
    ext_database_url_format = models.CharField(max_length=500)

    def __unicode__(self):
        """Returns the unicode representation of the external database."""
        return u'%s' % self.ext_database_name

class Curation_ExternalDatabase(models.Model):
    """'Through' table between Curation and ExternalDatabase tables."""

    # The curation object.
    curation = models.ForeignKey(Curation)

    # The external database object.
    external_database = models.ForeignKey(ExternalDatabase)

    # The external database accession number of the item.
    accession_number = models.CharField(max_length=500, null=False)

    def __unicode__(self):
        """Returns the unicode representation of the object."""
        return (u'curation: %d - xref: %s [%s]' %
                (self.curation.curation_id,
                 self.external_database.ext_database_name,
                 self.accession_number))

class NCBISubmission(models.Model):
    """NCBI submission table.

    The curated data in CollecTF is integrated into NCBI RefSeq database through
    periodic genome-specific submissions. This internal table keeps track of all
    the data that has been submitted to the NCBI RefSeq database.
    """
    # The time of the submission.
    submission_time = models.DateTimeField(auto_now_add=True)

    # The genome accession number that the submission is made for.
    genome_submitted_to = models.CharField(max_length=50)

    # The Curation_SiteInstance object that is submitted.
    curation_site_instance = models.ForeignKey('Curation_SiteInstance')

    class Meta:
        verbose_name = 'NCBI Submission'

