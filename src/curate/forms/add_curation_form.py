"""Form definitions for curation."""

import itertools
import re

from django import forms
from django.template.loader import get_template
from django.utils.safestring import mark_safe

from core import bioutils
from core import entrez_utils
from core.models import Curation
from core.models import ExperimentalTechnique
from core.models import ExternalDatabase
from core.models import Gene
from core.models import Genome
from core.models import GenomeSequence
from core.models import TF
from core.models import TFInstance
from core.models import Taxonomy
from curate import site_entry
from curate.forms import help_text


class PublicationForm(forms.Form):
    """Form for publication selection step.

    In this step, the curator is asked to select one of the papers assigned to
    him/her.
    """
    helptext = help_text.publication_form

    pub = forms.ChoiceField(
        widget=forms.RadioSelect(),
        label="Publications",
        help_text=helptext['pub'])

    no_data = forms.BooleanField(
        label="This paper contains no data.",
        required=False,
        help_text=helptext['no_data'])


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
        queryset=TF.objects.order_by('name'),
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
            Genome.objects.get(genome_accession=genome_accession)
        except Genome.DoesNotExist:
            # Get genome record
            try:
                record = entrez_utils.get_genome(genome_accession)
            except entrez_utils.EntrezException:
                raise forms.ValidationError("""
                Can not fetch genome record from NCBI. Check accession number.
                """)
            # Get taxonomy record
            try:
                entrez_utils.get_organism_taxon(record)
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
            new_genome(record, genes, species_taxon)

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
            TF_instance = TFInstance.objects.get(
                uniprot_accession=TF_accession)
        except TFInstance.DoesNotExist:
            # Create new TFInstance object in the database.
            try:
                entrez_utils.get_uniprot_TF(TF_accession)
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

            TFInstance.objects.create(
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
            genome = Genome.objects.get(genome_accession=genome_accession)
        except Genome.DoesNotExist:
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

        # Check if either TF_species or TF_species_same is filled
        if not (self.cleaned_data['TF_species'] or
                self.cleaned_data['TF_species_same']):
            self._errors['TF_species'] = self.error_class(
                ["Invalid TF species"])
        # Check if either site_species or site_species_same is filed
        if not (self.cleaned_data['site_species'] or
                self.cleaned_data['site_species_same']):
            self._errors['site_species'] = self.error_class(
                ["Invalid site species"])


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
        taxon, = filter(lambda x: x['Rank'] == rank, record['LineageEx'])
        lineage.append(taxon)
    parent = None
    for item in lineage:
        taxon_obj, _ = Taxonomy.objects.get_or_create(
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
    genome_seq = GenomeSequence.objects.create(
        sequence=str(genome_record.seq))
    # Create genome
    genome = Genome.objects.create(
        genome_accession=genome_record.id,
        genome_sequence=genome_seq,
        GC_content=bioutils.GC(genome_record.seq),
        gi=genome_record.annotations['gi'],
        organism=genome_record.annotations['organism'],
        chromosome=genome_record.features[0].qualifiers.get('chromosome', ''),
        taxonomy=species_taxon)
    # Create genes
    for gene in genes:
        Gene.objects.create(genome=genome, **gene)
    return genome


class TechniquesForm(forms.Form):
    """Form to enter experimental techniques used to identify TFBS.

    Curators are asked to specify ALL techniques that are used to identify the
    sites that they plan to report in the curation. The total number of
    different techniques will be used to create and populate the columns in the
    final form.

    Curators are also asked to provide a brief description of the experimental
    setup for the sites reported. In this new setup, the description should
    comprise the setup used for all the sites that are going to be reported,
    even though they might used different techniques. For instance: "The
    binding motif for XX was identified through phylogenetic footprinting and
    directed-mutagenesis + EMSA on the promoter of gene YY. Researchers then
    performed a computer search of sites with one mismatch in the genome. Of
    those identified, they verified X through EMSA. Regulatory activity for
    four of the sites was assessed with beta-gal assays.". In brief, the
    curator is expected to provide a concise and logical summary of the
    experimental process leading to the identification of all reported sites
    and the determination (if any) of their regulatory activity.

    Curators, as before, will be prompted to specify any external DBs where
    high-throughput data might be stored (e.g. array data on GEO).

    Curators will also be given the option to specify whether the TF is shown
    to interact with another protein/ligand that influences binding (and
    optionally add notes on that [pop-up]).
    """

    NUM_EXTRA_DB_FIELDS = 5

    def __init__(self, *args, **kwargs):
        """Overrides initialization"""
        super(TechniquesForm, self).__init__(*args, **kwargs)
        help_dict = help_text.techniques_form
        # Extra external-database fields
        external_db_type_choices = [(None, 'None')]
        for db in ExternalDatabase.objects.all():
            external_db_type_choices.append(
                (db.ext_database_id, db.ext_database_name))
        # Create extra fields
        for i in range(self.NUM_EXTRA_DB_FIELDS):
            self.fields['external_db_type_%d' % i] = forms.ChoiceField(
                choices=external_db_type_choices,
                required=False,
                label="External DB type [%d]" % (i+1),
                help_text=help_dict['external_db_type'])
            self.fields['external_db_accession_%d' % i] = forms.CharField(
                required=False,
                label="External DB accession number [%d]" % (i+1),
                help_text=help_dict['external_db_accession'])

    help_dict = help_text.techniques_form

    # Generate techniques field by getting available techniques from db
    template = get_template('experimental_technique_field.html')
    choices = [(technique.technique_id,
                template.render({'technique': technique}))
               for technique in ExperimentalTechnique.objects.order_by('name')]

    techniques = forms.MultipleChoiceField(
        choices=choices, label="Techniques", help_text=help_dict['techniques'],
        widget=forms.CheckboxSelectMultiple())

    experimental_process = forms.CharField(
        widget=forms.Textarea,
        label="Experimental process",
        help_text=help_dict['experimental_process'])

    # Does TF interact with any other protein/ligand that influences binding?
    forms_complex = forms.BooleanField(
        required=False,
        label="""
        The manuscript reports that TF forms complex with other proteins for
        binding with reported sites.""")

    complex_notes = forms.CharField(
        widget=forms.Textarea, required=False, label="Notes",
        help_text=help_dict['complex_notes'])

    # External database links
    has_external_db = forms.BooleanField(
        required=False,
        label="""The manuscript reports high-throughput data from an external
        database. (You can report up to 5 external resources.)""")


class SiteEntryForm(forms.Form):
    """Form for reporting sites.

    Curators are first asked to specify whether the sites are: motif associated
    [MA], variable motif associated (e.g. variable spacing,
    inverting... anything that is not gapless alignment) [VMA] or non-motif
    associated (no specific sequence pattern has been determined as the binding
    motif) [NMA].

    Curators are later asked to enter the sites. They can do this in two major
    formats: (1) Sequence-based (e.g. CTGTTGCACGT) (2) Coordinate-based
    (e.g. 12312 12323)

    If the user has not checked "This is the same strand used in the paper" in
    the first page, the system will ask the user to verify that the coordinates
    they enter refer to the NCBI strand.

    - Curators can also add quantitative data to either format.
    - All fields must be separated by either space or tab
    - The system will recognize the entry format and the presence of
    quantitative data once Next is clicked. If there is quantitative data, the
    system will prompt the user for a brief description of the field.
    """

    help_dict = help_text.site_entry_form
    # Type of sites to be entered
    # The curator is able to choose one of the available motifs or create a new
    # one. Sites can also be curated as either variable-motif-associated or
    # non-motif-associated.
    # To be populated dynamically
    site_type = forms.ChoiceField(
        choices=(),
        widget=forms.RadioSelect,
        required=True,
        label="Site type",
        help_text=mark_safe(help_dict['site_type']))

    sites = forms.CharField(
        required=True, widget=forms.Textarea,
        label="Sites",
        help_text=mark_safe(help_dict['sites']))

    quantitative_data_format = forms.CharField(
        required=False,
        label="Quantitative data format",
        help_text=help_dict['quantitative_data_format'])

    # Following fields will be visible only if the submission is
    # high-throughput
    peaks = forms.CharField(
        widget=forms.Textarea,
        label="High-throughput sequences",
        help_text=help_dict['peaks'])

    assay_conditions = forms.CharField(
        label="Assay conditions",
        help_text=help_dict['assay_conditions'],
        widget=forms.Textarea)

    method_notes = forms.CharField(
        label="Method notes",
        help_text=help_dict['method_notes'],
        widget=forms.Textarea)

    peak_techniques = forms.MultipleChoiceField(
        label="Techniques used to identify high-throughput data",
        help_text=help_dict['peak_techniques'],
        choices=(),
        widget=forms.CheckboxSelectMultiple)

    def verify_coordinates(self, coor_a, coor_b):
        """Verifies if coordinates are positive integers."""
        try:
            ca, cb = int(coor_a), int(coor_b)
            return ca, cb if ca <= 0 or cb <= 0 else None
        except ValueError:
            return None

    def verify_float(self, s):
        """Checks if the string is in float format."""
        try:
            float(s)
            return True
        except ValueError:
            return None

    def verify_only_sites(self, sites_cd):
        """Checks if all site sequences are valid."""
        try:
            # Check if it is FASTA format.
            if sites_cd.startswith('>'):
                site_entry.parse_fasta(sites_cd)
            else:
                lines = [line.split() for line in sites_cd.split('\n')]
                site_entry.parse_seq(sites_cd)
        except:
            raise forms.ValidationError("Ambiguous DNA sequence.")
        return sites_cd

    def verify_only_coordinates(self, sites_cd):
        """Checks if all coordinates are valid."""
        coordinates = [re.split('[\t ]+', line)
                       for line in re.split('[\r\n]+', sites_cd)]

        msg = "Invalid coordinate format."
        for instance in coordinates:
            if len(instance) != 2:
                raise forms.ValidationError(msg)
            if not self.verify_coordinates(instance[0], instance[1]):
                raise forms.ValidationError(msg)
        return sites_cd

    def verify_only_sites_and_values(self, sites_cd):
        """Verifies entries with a sequence and quantitative value per site."""
        lines = [line.split() for line in sites_cd.split('\n')]
        if (all(len(line) == 2 for line in lines) and
            all(nuc in 'ACTG' for line in lines for nuc in line[0]) and
            all(self.verify_float(line[1]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("Invalid format")

        return sites_cd

    def verify_only_coordinates_and_values(self, sites_cd):
        """Verifies site entries with coordinates and quantitative values."""
        lines = [line.split() for line in sites_cd.split('\n')]
        if (all(len(line) == 3 for line in lines) and
            all(self.verify_coordinates(line[0], line[1]) for line in lines) and
            all(self.verify_float(line[2]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("Invalid input format.")
        return sites_cd

    def check_quantitative_data_format(self, qval_data_format):
        """Checks if the quantitative-data-format is filled."""
        if not qval_data_format:
            msg = """
            The reported sites have associated quantitative values. Please enter
            the quantitative data format."""
            self._errors["__all__"] = self.error_class([msg])

    def clean_sites(self):
        """Validates sites field."""
        cd = self.cleaned_data['sites'].strip().upper()
        lines = [re.split('[\t ]+', line.strip())
                 for line in re.split('[\r\n]+', cd.strip())]
        sites_cd = '\n'.join(' '.join(wds for wds in line) for line in lines)

        # By default, don't check the qval data format.
        self.check_qval_data_format = False
        # Check if it is in FASTA format.
        if sites_cd[0].startswith('>'):
            return self.verify_only_sites(sites_cd)
        if len(lines[0]) == 1:
            return self.verify_only_sites(sites_cd)
        elif len(lines[0]) == 2:
            if lines[0][0].isalpha():
                self.check_qval_data_format = True
                return self.verify_only_sites_and_values(sites_cd)
            else:
                return self.verify_only_coordinates(sites_cd)
        elif len(lines[0]) == 3:
            self.check_qval_data_format = True
            return self.verify_only_coordinates_and_values(sites_cd)

        # Otherwise return validation error.
        raise forms.ValidationError("Invalid format.")

    def clean_peaks(self):
        """Validates peaks field."""
        cd = self.cleaned_data['peaks'].strip().upper()
        lines = [re.split('[\t ]+', line.strip())
                 for line in re.split('[\r\n]+', cd.strip())]
        peaks_cd = '\n'.join(' '.join(wds for wds in line) for line in lines)

        # By default, don't check the qval data format.
        self.check_qval_data_format = False
        if len(lines[0]) == 1:  # sequence format
            # Same process with site verification.
            return self.verify_only_sites(peaks_cd)
        elif len(lines[0]) == 2:
            # Either coordinate format or sequence with qvals
            if lines[0][0].isalpha():  # seq with qvals
                self.check_qval_data_format = True
                return self.verify_only_sites_and_values(peaks_cd)
            else:
                return self.verify_only_coordinates(peaks_cd)
        elif len(lines[0]) == 3:  # coordinates with qval
            self.check_qval_data_format = True
            return self.verify_only_coordinates_and_values(peaks_cd)

        raise forms.ValidationError("Invalid format")

    def clean(self):
        """Cleans form fields."""
        # If the site field contains quantitative data, make sure the format
        # field is filled.
        if 'sites' in self.cleaned_data and self.check_qval_data_format:
            self.check_quantitative_data_format(
                self.cleaned_data.get('quantitative_data_format', None))
        return self.cleaned_data


class SiteExactMatchForm(forms.Form):
    """Form to select and match reported sites in the genome.

    This form displays only exact matches (i.e. ones that is present in genome
    exactly). If the site is not found in the genome 'exactly', it is searched
    'softly' and presented in the next form, SiteSoftMatchForm.
    """
    # All form fields are created dynamically. No static field def here
    def clean(self):
        """Cleans fields."""
        return self.cleaned_data


class SiteSoftMatchForm(forms.Form):
    """Form displaying results of 'soft' search.

    Match sites are not exactly same with the query sequence, but similar."""
    # No static def here either.
    pass


class SiteAnnotationForm(forms.Form):
    """Form asking the curator to fill the information regarding each site.

    In particular, the user can:
    - visualize again site information, including chromosome,
    - toggle the graphical view (off by default)
    - edit the qualitative values
    - specify which techniques were used to determine this site
    - define the site as repressed, activated or whether the effect of TF on
      the site is unknown (not-determined)

    Each site and its set of fields that can be edited are represented as a
    form (one form per site instance). The abstraction to work multiple forms
    in one page is achieved via Django FormSets
    (https://docs.djangoproject.com/en/1.6/topics/forms/formsets/).
    """
    def clean(self):
        """Cleans form fields."""
        # TODO(sefa): Check if at least one experimental techniques is selected
        # for each site.
        return self.cleaned_data


class GeneRegulationForm(forms.Form):
    """Form to enter gene regulation data.

    This form is displayed after SiteSoftMatchForm. After the curator selects
    site equivalent for each reported site, in this form, surrounding genes to
    the site are displayed. For each gene, it can be (un)checked whether site
    regualates gene (or not).

    Like the previous two forms (SiteExactMatchForm and SiteSoftMatchForm), all
    fields in this form are created dynamically, based on which genome
    positions are selected in the previous two forms as site equivalents.
    """
    pass


class CurationReviewForm(forms.Form):
    """Form to review all the data entered so far.

    This is the last step before submitting curation.
    """
    help_dict = help_text.curation_review_form

    choices = ((None, 'None'),) + Curation.REVISION_REASONS
    revision_reasons = forms.ChoiceField(
        choices=choices,
        label="Revision required",
        help_text=help_dict['revision_reasons'],
        required=False)

    confidence = forms.BooleanField(
        required=False,
        label="I am confident of the results reported in this manuscript.",
        help_text=help_dict['confidence'])

    paper_complete = forms.BooleanField(
        required=False,
        label="Curation for this paper is complete.",
        help_text=help_dict['paper_complete'])

    notes = forms.CharField(widget=forms.Textarea,
                            required=False,
                            label="Notes",
                            help_text=help_dict['notes'])

    confirm = forms.BooleanField(required=True,
                                 label="I want to submit this curation",
                                 help_text=help_dict['confirm'])
