import itertools

from django import forms

from core import models


class PublicationForm(forms.Form):
    """Publication selection form.

    The curators are asked to select one of the papers assigned to them.
    """
    publication = forms.ChoiceField(
        label="Publications",
        widget=forms.RadioSelect())

    no_data = forms.BooleanField(
        label="This paper contains no data.",
        required=False,
        help_text="""
        Check this button if, after examining the paper, you find that the
        paper does not have data on binding sites. Checking this button will
        mark the paper as having no binding site data and set it to the
        'curation complete' status. Also, the curation process will be ended as
        the paper has no data to be curated.""")


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
                required=False,
                widget=forms.HiddenInput())

        # Order fields
        self.order_fields(
            ['TF',
             'genome_accession'] +
            ['genome_accession_%d' % i
             for i in range(self.NUM_EXTRA_GENOME_FIELDS)] +
            ['site_species_same',
             'TF_accession',
             'TF_refseq_accession'] +
            list(itertools.chain(
                ['TF_accession_%d' % i
                 for i in range(self.NUM_EXTRA_TF_FIELDS)],
                ['TF_refseq_accession_%d' % i
                 for i in range(self.NUM_EXTRA_TF_FIELDS)])) +
            ['TF_species_same',
             'site_species',
             'TF_species',
             'contains_promoter_data',
             'contains_expression_data'])

    TF = forms.ModelChoiceField(
        queryset=models.TF.objects.order_by('name'),
        label="TF",
        help_text="""
        Select the transcription factor you are curating on from list. If not
        in list, please add the corresponding TF/TF-family using Data
        submission menu.""")

    genome_accession = forms.CharField(
        label="Genome NCBI accession number",
        help_text="""
        Paste the NCBI GenBank genome accession number for the species closest
        to the reported species/strain.  (e.g. <code>NC_000913.2</code>) You
        can add more than one chromosome.""")

    site_species_same = forms.BooleanField(
        required=False,
        label="""
        This is the exact same strain as reported in the manuscript for the
        sites.""",
        help_text="""
        Check if the reported strain and selected RefSeq strain are same.""")

    TF_accession = forms.CharField(
        label="TF UniProt accession number",
        help_text="""
        Enter the UniProt TF protein accession number for the species closest
        to the reported species/strain.  (e.g. <code>Q6GZX3</code>) You can add
        more than one TF.""")


    TF_refseq_accession = forms.CharField(
        label="TF RefSeq accession number",
        help_text="""
        Enter the RefS eq TF protein accession number for the species closest to
        the reported species/strain.  (e.g. <code>NP_799324</code>) You can add
        more than one TF.""",
        widget=forms.HiddenInput)

    TF_species_same = forms.BooleanField(
        required=False,
        label="""
        This is the exact same strain as reported in the manuscript for the
        TF.""",
        help_text="""
        Check if the reported strain and selected RefSeq strain are same for
        the TF.  If it is checked and TF accession has <code>WP_</code> prefix,
        organism of origin for reported TF will be set to organism of entered
        genome.""")

    site_species = forms.CharField(
        label="Organism TF binding sites are reported in",
        required=False,
        help_text="""
        If the work you are reporting uses a strain different from the selected
        RefSeq genome, please type/paste the original strain
        (e.g. <code>Pseudomonas putida plasmid pEST1226</code>). This allows us
        to keep track of the correspondence between reported and mapped
        strains.""")

    TF_species = forms.CharField(
        label="Organism of origin for reported TF",
        required=False,
        help_text="""
        If the work you are reporting uses a strain different from the selected
        RefSeq genome, please type/paste the original strain
        (e.g. <code>Pseudomonas sp. ADP</code>). This allows us to keep track
        of the correspondence between reported and mapped strains.""")

    # Publication-related fields
    contains_promoter_data = forms.BooleanField(
        required=False,
        label="The manuscript contains promoter information",
        help_text="""
        Check if the paper provides experimental data on the structure and
        sequence of a TF-regulated promoter""")

    contains_expression_data = forms.BooleanField(
        required=False,
        label="The manuscript contains expression data",
        help_text="""
        Check if the paper provides experimental support for TF-mediated
        regulation of genes.  Please make sure that this field is checked if
        you plan to report differential gene expression associated with TF
        activity.""")




    def clean_genome_accession(self):
        """Cleans the genome accession field."""
