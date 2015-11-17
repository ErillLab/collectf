import itertools

from django import forms

from core import models
from core import entrez_utils
from core import bioutils
from . import help_text


class GenomeForm(forms.Form):
    """Genome and TF accession number form."""
    NUM_EXTRA_GENOME_FIELDS = 3
    NUM_EXTRA_TF_FIELDS = 3

    def __init__(self, *args, **kwargs):
        super(GenomeForm, self).__init__(*args, **kwargs)

        # Extra genome accession fields
        for i in range(self.NUM_EXTRA_GENOME_FIELDS):
            self.fields['genome_accession_%d' % i] = forms.CharField(
                label="Genome NCBI accession number (%d)" % i,
                required=False)

        # Extra TF accession fields
        self.num_extra_TF_fields = 3
        for i in range(self.NUM_EXTRA_TF_FIELDS):
            self.fields['TF_accession_%d' % i] = forms.CharField(
                label="TF UniProt accession number (%d)" % i,
                required=False)
            self.fields['TF_refseq_accession_%d' % i] = forms.CharField(
                label="TF RefSeq accession number (%d)" % i,
                required=False)

        # Order fields
        self.order_fields(
            ['TF',
             'genome_accession'] +
            ['genome_accession_%d' % i
             for i in range(self.NUM_EXTRA_GENOME_FIELDS)] +
            ['site_species_same',
             'TF_accession',
             'TF_refseq_accession'] +
            list(itertools.chain(*zip(
                ['TF_accession_%d' % i
                 for i in range(self.NUM_EXTRA_TF_FIELDS)],
                ['TF_refseq_accession_%d' % i
                 for i in range(self.NUM_EXTRA_TF_FIELDS)]))) +
            ['TF_species_same',
             'site_species',
             'TF_species',
             'contains_promoter_data',
             'contains_expression_data'])

    help_text_dict = help_text.genome_form

    TF = forms.ModelChoiceField(
        queryset=models.TF.objects.order_by('name'),
        label="TF",
        help_text=help_text_dict['TF'])

    genome_accession = forms.CharField(
        label="Genome NCBI accession number",
        help_text=help_text_dict['genome_accession'])

    site_species_same = forms.BooleanField(
        required=False,
        label="""
        This is the exact same strain as reported in the manuscript for the
        sites.""",
        help_text=help_text_dict['site_species_same'])

    TF_accession = forms.CharField(
        label="TF UniProt accession number",
        help_text=help_text_dict['TF_accession'])

    TF_refseq_accession = forms.CharField(
        label="TF RefSeq accession number",
        help_text=help_text_dict['TF_refseq_accession'])

    TF_species_same = forms.BooleanField(
        required=False,
        label="""
        This is the exact same strain as reported in the manuscript for the
        TF.""",
        help_text=help_text_dict['TF_species_same'])

    site_species = forms.CharField(
        label="Organism TF binding sites are reported in",
        required=False,
        help_text=help_text_dict['site_species'])

    TF_species = forms.CharField(
        label="Organism of origin for reported TF",
        required=False,
        help_text=help_text_dict['TF_species'])

    # Publication-related fields
    contains_promoter_data = forms.BooleanField(
        required=False,
        label="The manuscript contains promoter information",
        help_text=help_text_dict['contains_promoter_data'])

    contains_expression_data = forms.BooleanField(
        required=False,
        label="The manuscript contains expression data",
        help_text=help_text_dict['contains_expression_data'])

    def clean_genome_accession_helper(self, genome_accession):
        """Checks if the genome accession field is valid.

        If the genome accession is valid and not in the database, downloads
        genome sequences and list of genes from NCBI database and adds them to
        the CollecTF database.

        Gets called for all genome accesion fields.
        """
        genome_accession = genome_accession.strip()
        if '.' not in genome_accession:
            raise forms.ValidationError(
                "Please enter RefSeq accession number with the version number.")
        if not genome_accession.startswith('NC_'):
            raise forms.ValidationError(
                "RefSeq genome accession number should start with 'NC_'")

        try:
            models.Genome.objects.get(genome_accession=genome_accession)
        except models.Genome.DoesNotExist:
            # Get genome record
            try:
                record = entrez_utils.get_genome(genome_accession)
            except entrez_utils.EntrezException:
                raise forms.ValidationError("""
                Can not fetch genome record from NCBI. Check accession number.
                """)
            # Get taxonomy record
            try:
                strain_tax = entrez_utils.get_organism_taxon(record)
            except entrez_utils.EntrezException:
                raise forms.ValidationError(
                    "Can not fetch strain taxonomy information.")

            # Get genes
            try:
                genes = entrez_utils.get_genes(record)
            except entrez_utils.EntrezException:
                raise forms.ValidationError("""
                Can't retrieve list of genes. Check genome accession
                number.""")

            # Create genome object and genes.
            species_taxon = new_taxonomy(record)
            genome = new_genome(record, genes, species_taxon)
            
        return genome_accession

    def clean_genome_accession(self):
        """Cleans genome accession field."""
        genome_accession = self.cleaned_data['genome_accession']
        return self.clean_genome_accession_helper(genome_accession)

    def check_TF_accession(self, TF_accession, TF_refseq_accession):
        """Checks if the entered TF accession number is valid."""
        TF = self.cleaned_data.get('TF')
        if not TF:
            return
        try:
            TF_instance = models.TFInstance.objects.get(
                uniprot_accession=TF_accession)
        except models.TFInstance.DoesNotExist:
            # Create new TFInstance object in the database.
            try:
                TF_record = entrez_utils.get_uniprot_TF(TF_accession)
            except entrez_utils.EntrezException:
                raise forms.ValidationError("""
                Can not fetch protein record from UniProt. Check accession
                numbers.""")

            if not TF_refseq_accession:
                raise forms.ValidationError(
                    "Enter TF RefSeq accession number.")

            # Check if entered/retrieved WP accession is valid
            TF_refseq_accession = TF_refseq_accession.strip().split('.')[0]
            if not TF_refseq_accession.startswith('WP_'):
                raise forms.ValidationError("Invalid RefSeq accession.")
            try:
                refseq_record = entrez_utils.get_refseq_TF(TF_refseq_accession)
            except entrez_utils.EntrezException:
                raise forms.ValidationError("Invalid RefSeq accession.")
                
            models.TFInstance.objects.create(
                uniprot_accession=TF_accession,
                refseq_accession=TF_refseq_accession,
                TF=TF,
                description=refseq_record.description)
        else:
            # Check if selected TF matches with the TF_instance's TF.
            if TF != TF_instance.TF:
                raise forms.ValidationError(
                    "It seems that %s is a %s but you selected %s." %
                    (TF_accession, TF_instance.TF.name, TF.name))

    def clean_TF_accession_helper(self, TF_accession):
        return TF_accession.strip()
    
    def clean_TF_accession(self):
        return self.clean_TF_accession_helper(
            self.cleaned_data['TF_accession'])

    def clean_TF_refseq_accession_helper(self, TF_refseq_accession):
        refseq_accession = TF_refseq_accession.strip().split('.')[0]
        if refseq_accession.startswith('WP'):
            try:
                entrez_utils.get_refseq_TF(refseq_accession)
                return refseq_accession
            except entrez_utils.EntrezException:
                pass
        raise forms.ValidationError("Invalid RefSeq accession.")

    def clean_TF_refseq_accession(self):
        return self.clean_TF_refseq_accession_helper(
            self.cleaned_data['TF_refseq_accession'])

    def clean_species(self, field):
        """Helper function for clean_TF_species and clean_site_species.
        
        When TF_species_same or site_species_same fields are selected, returns
        the organism information from the entered genome accession number.
        """
        genome_accession = self.cleaned_data.get('genome_accession')
        if not genome_accession:
            return ''
        try:
            genome = models.Genome.objects.get(
                genome_accession=genome_accession)
        except models.Genome.DoesNotExist:
            raise forms.ValidationError("Invalid genome accession number")
        return genome.organism

    def clean_TF_species(self):
        """Cleans TF_species field.
        
        If TF_species_same field is selected, assigns species data to the
        cleaned data. In that case, it is important that clean_genome is called
        BEFORE clean_TF_species, because the genome is needed to extract
        species information.
        """
        return (self.clean_species('TF_species')
                if self.cleaned_data['TF_species_same']
                else self.cleaned_data['TF_species'])

    def clean_site_species(self):
        """Cleans site_species field.
        
        If site_species_same field is selected, assigns species data to the
        cleaned data.
        """
        return (self.clean_species('site_species')
                if self.cleaned_data['site_species_same']
                else self.cleaned_data['site_species'])

    def clean(self):
        # Clean extra genome accession fields
        for i in range(self.NUM_EXTRA_GENOME_FIELDS):
            field = 'genome_accession_%d' % i
            if self.cleaned_data[field]:
                self.cleaned_data[field] = self.clean_genome_accession_helper(
                    self.cleaned_data[field])

        # Clean all TF accesion fields
        if ('TF_accession' in self.cleaned_data and
            'TF_refseq_accession' in self.cleaned_data):
            self.check_TF_accession(
                self.cleaned_data['TF_accession'],
                self.cleaned_data['TF_refseq_accession'])
        for i in range(self.NUM_EXTRA_TF_FIELDS):
            uniprot = 'TF_accession_%d' % i
            refseq = 'TF_refseq_accession_%d' % i
            if not (self.cleaned_data[uniprot] and self.cleaned_data[refseq]):
                continue
            self.cleaned_data[uniprot] = self.clean_TF_accession_helper(
                self.cleaned_data[uniprot])
            self.cleaned_data[refseq] = self.clean_TF_refseq_accession_helper(
                self.cleaned_data[refseq])
            self.check_TF_accession(self.cleaned_data[uniprot],
                                    self.cleaned_data[refseq])


def new_taxonomy(genome_record):
    """Creates all taxonomy items up to Bacteria, given the genome record.

    Returns the Taxonomy object at the species level.
    """
    record = entrez_utils.get_taxonomy(genome_record)
    # Sort lineage
    # The lineage that NCBI returns seems already sorted, but to be safe, sort
    # it here.
    lineage = []
    for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
        taxon, = filter(lambda x: x['Rank']==rank, record['LineageEx'])
        lineage.append(taxon)
    parent = None
    for item in lineage:
        taxon_obj, _ = models.Taxonomy.objects.get_or_create(
            rank=item['Rank'],
            taxonomy_id=item['TaxId'],
            name=item['ScientificName'],
            parent=parent)
        parent = taxon_obj
        
    return taxon_obj

        
def new_genome(genome_record, genes, species_taxon):
    """Creates a models.Genome object from the given genome record.

    Returns the created genome object."""
    # Create genome sequence
    genome_seq = models.GenomeSequence.objects.create(
        sequence = str(genome_record.seq))
    # Create genome
    genome = models.Genome.objects.create(
        genome_accession=genome_record.id,
        genome_sequence=genome_seq,
        GC_content=bioutils.GC(genome_record.seq),
        gi=genome_record.annotations['gi'],
        organism=genome_record.annotations['organism'],
        chromosome=genome_record.features[0].qualifiers.get('chromosome', ''),
        taxonomy=species_taxon)
    # Create genes
    for gene in genes:
        models.Gene.objects.create(genome=genome, **gene)
    return genome

    
