from django import forms
from django.forms.formsets import BaseFormSet
from base import bioutils
# import all models
from base.models import *
from curate import create_object
from curate import site_entry
import re
import help_texts
from django.utils.safestring import mark_safe
from collectf import settings
from django.template.loader import get_template
from django.template import Context


class PublicationForm(forms.Form):
    """Form for publication selection step.  In this step, the user is asked to
    select one of the papers assigned to him/her."""
    ht = help_texts.publication_form
    pub = forms.ChoiceField(widget=forms.RadioSelect(),
                            label="Publications", help_text=ht['pub'])
    no_data = forms.BooleanField(label="This paper contains no data.",
                                 required=False, help_text=ht['no_data'])


class GenomeForm(forms.Form):
    """Form for submission of genome and TF accession numbers and others.

    - The user is asked for the TF and and its reported structure (if multiple
      structures are reported in a paper, this is one of the few things that
      will still require multiple submissions, since motifs are likely to
      change).
    - Then the user is asked to enter at least one genome and TF accession
      number. The system checks for all genomes and for all TFs, but separately,
      that they come from the same taxonomy ID. This allows the user to deal
      automagically with bacteria with multiple chromosomes, and it allows the
      user to define multi-unit TFs, such as heterodimers like IHF.
    - The user is finally asked for paper reported site and TF species.
    - The user is also allowed to edit the two manuscript fields (optional) that
      specify whether the paper contains expression data and promoter
      information.
      """

    def __init__(self, *args, **kwargs):
        """Override initialization"""
        super(GenomeForm, self).__init__(*args, **kwargs)
        num_genome_fields = settings.NUMBER_OF_GENOME_ACCESSION_FIELDS
        num_TF_fields = settings.NUMBER_OF_TF_ACCESSION_FIELDS
        # Extra genome accession fields
        for i in xrange(1, num_genome_fields):
            self.fields['genome_accession_%d' % i] = forms.CharField(
                label="Genome NCBI accession number [%d]" % i,
                required=False)
        # Extra TF accession fields
        for i in xrange(1, num_TF_fields):
            self.fields['TF_accession_%d' % i] = forms.CharField(
                label="TF accession number [%d]" % i,
                required=False)
        # Change the order of fields
        current_order = self.fields.keyOrder
        self.fields.keyOrder = (['TF', 'genome_accession'] +
                                ['genome_accession_%d' % i
                                 for i in xrange(1, num_genome_fields)] +
                                ['site_species_same'] +
                                ['TF_accession'] +
                                [('TF_accession_%d' % i)
                                 for i in xrange(1, num_TF_fields)] +
                                ['TF_species_same',
                                 'site_species',
                                 'TF_species'] +
                                ['contains_promoter_data',
                                 'contains_expression_data'])

    help_dict = help_texts.genome_form
    # TF field
    TF = forms.ModelChoiceField(queryset=TF.objects.order_by('name'),
                                label="TF",
                                help_text=help_dict['TF'])

    # Genome accession number(s)
    # The last two genome accession fields will be hidden by default, but be
    # able to shown by the curator for data entry.
    genome_accession = forms.CharField(label="Genome NCBI accession number",
                                       help_text=help_dict['genome_accession'])

    # Checked if site species is the same with the reported genome
    site_species_same = forms.BooleanField(
        required=False,
        label="""This is the exact same strain as reported in the manuscript for
        the sites.""",
        help_text=help_dict['site_species_same'])
    # TF accession number
    TF_accession = forms.CharField(label="TF accession number",
                                   help_text=help_dict['TF_accession'])

    # Checked if TF species is the same with the reported genome
    TF_species_same = forms.BooleanField(
        required=False,
        label="""This is the exact same strain as
        reported in the manuscript for the TF.""",
        help_text=help_dict['TF_species_same'])

    site_species = forms.CharField(
        label="Organism TF binding sites are reported in",
        required=False,
        help_text=help_dict['site_species'])

    TF_species = forms.CharField(label="Organism of origin for reported TF",
                                 required=False,
                                 help_text=help_dict['TF_species'])

    # manuscript fields edit.
    contains_promoter_data = forms.BooleanField(
        required=False,
        label="The manuscript contains promoter information",
        help_text=help_dict['contains_promoter_data'])

    contains_expression_data = forms.BooleanField(
        required=False,
        label="The manuscript contains expression data",
        help_text=help_dict['contains_expression_data'])

    def clean_genome_accession_helper(self, genome_accession):
        """Check if the entered genome accession number is valid"""
        if '.' not in genome_accession:
            msg = """Please enter RefSeq accession number with the version
            number."""
            raise forms.ValidationError(msg)
        if not genome_accession.startswith('NC_'):
            msg = "RefSeq genome accession number should start with 'NC_'"
            raise forms.ValidationError(msg)

        try: # to retrieve genome from database
            g = Genome.objects.get(genome_accession=genome_accession)
        except Genome.DoesNotExist: # try to retrieve it from NCBI database
            genome_record = bioutils.get_genome(genome_accession)
            strain_tax = bioutils.get_organism_taxon(genome_record)
            if not genome_record:
                msg = """Can not fetch genome record from NCBI.
                Check accession number."""
                raise forms.ValidationError(msg)
            if not strain_tax:
                msg = "Can not fetch strain taxonomy information."
                raise forms.ValidationError(msg)
            # At this point, things should be fine, the next step is to create
            # genome object and genes.
            if not create_object.make_genome(genome_record, strain_tax):
                msg = """Can't retrieve genome or list of genes.
                Does the genome have a RefSeq accession number?"""
                raise forms.ValidationError(msg)

        return genome_accession

    def clean_genome_accession(self):
        """Clean genome accession fields. Check if they are valid. If they are
        valid and not available in the database, download genome sequences and
        list of genes from NCBI database and add them to the CollecTF
        database. Do the validity check for all genome accession fields (if
        there are more than one)."""
        genome_accession = self.cleaned_data['genome_accession'].strip()
        return self.clean_genome_accession_helper(genome_accession)

    def clean_TF_accession_helper(self, TF_accession):
        """Check if the entered TF accession number is valid"""
        try:
            TF_accession = TF_accession.split('.')[0]
            if not (TF_accession.startswith('NP_') or
                    TF_accession.startswith('YP_') or
                    TF_accession.startswith('WP_')):
                msg = """TF accession number must start with 'NP_', 'YP_' or 'WP_'."""
                raise forms.ValidationError(msg)
            TF_instance = TFInstance.objects.get(protein_accession=TF_accession)
        except TFInstance.DoesNotExist:
            TF_record = bioutils.get_TF(TF_accession)
            if not TF_record:
                msg = """Can not fetch protein record from NCBI. Check accession
                number."""
                raise forms.ValidationError(msg)
            # Create TF instance object
            create_object.make_TF_instance(TF_record)
        return TF_accession

    def clean_TF_accession(self):
        """Clean TF accession fields. Check if they are valid. If they are valid
        and not available in the database, download them from the NCBI database
        and add to the CollecTF . Do the validity check for all TF accession
        fields (if there are more than one)."""
        TF_accession = self.cleaned_data['TF_accession'].strip()

        return self.clean_TF_accession_helper(TF_accession)

    def clean_species(self, field):
        """Helper function for clean_TF_species and clean_site_species. When
        TF_species_same or site_species_same fields are selected, it returns the
        organism information from the entered genome accession number."""
        if 'genome_accession' not in self.cleaned_data:
            return
        genome_accession = self.cleaned_data['genome_accession']
        genome = Genome.objects.get(genome_accession=genome_accession)
        if not genome:
            msg = "Invalid genome accession number"
            self._errors[field] = self.error_class([msg])
        return genome.organism

    def clean_TF_species(self):
        """Clean TF_species field. If TF_species_same field is selected, assign
        species data to the cleaned data. In that case, it is important that
        clean_genome is called BEFORE clean_TF_species, because the genome is
        needed to extract species information."""
        return (self.clean_species('TF_species')
                if self.cleaned_data['TF_species_same']
                else self.cleaned_data['TF_species'])

    def clean_site_species(self):
        """Clean site_species field. If site_species_same field is selected,
        assign species data to the cleaned data."""
        return (self.clean_species('site_species')
                if self.cleaned_data['site_species_same']
                else self.cleaned_data['site_species'])

    def check_genome_accession_origin(self):
        """Check if all genome accession numbers belong to the same taxonomy
        ID"""
        cd = self.cleaned_data
        if 'genome_accession' not in cd:
            return
        # Check extra genome-accession fields
        genome_accessions = [cd['genome_accession']]
        for i in xrange(1, settings.NUMBER_OF_GENOME_ACCESSION_FIELDS):
            field_name = 'genome_accession_%d' % i
            if cd.get(field_name, None):
                self.clean_genome_accession_helper(cd[field_name].strip())
                genome_accessions.append(cd[field_name].strip())
        # Get all genomes from the database.
        genomes = [Genome.objects.get(genome_accession=acc)
                   for acc in genome_accessions]
        all_same = lambda items: all(x == items[0] for x in items)
        if not all_same([genome.organism for genome in genomes]):
            msg = "Genome accession numbers are not from the same taxonomy ID."
            self._errors['genome_accession'] = self.error_class([msg])

    def check_TF_accession_origin(self):
        """Check if all TF accession fields belong to the same taxonomy ID"""
        try:
            cd = self.cleaned_data
            if 'TF_accession' not in cd:
                return
            TF_accessions = [cd['TF_accession']]
            if len(TF_accessions) > 1:
                for i in xrange(settings.NUMBER_OF_TF_ACCESSION_FIELDS):
                    field_name = 'TF_accession_%d' % i
                    if cd.get(field_name, None):
                        self.clean_TF_accession_helper(cd[field_name].strip())
                        TF_accessions.append(cd[field_name].strip())
                # Check if all TF accession numbers come from the same organism
                all_same = lambda items: all(x == items[0] for x in items)
                if not all_same([bioutils.TF_accession_to_org_taxon(acc)
                                 for acc in TF_accessions]):
                    msg = "TF accession numbers are not from the same taxonomy ID."
                    self._errors['TF_accession'] = self.error_class([msg])

        except:
            msg = """Failed to validate TF accession numbers (can not fetch
            records from NCBI)"""
            self._errors['TF_accession'] = self.error_class([msg])

    def clean(self):
        """All other clean operations"""
        # Check if either TF_species or TF_species_same is filled
        cd = self.cleaned_data
        if not (cd['TF_species'] or cd['TF_species_same']):
            self._errors['TF_species'] = self.error_class(["Invalid TF species"])
        # Check if either site_species or site_species_same is filed
        if not (cd['site_species'] or cd['site_species_same']):
            self._errors['site_species'] = self.error_class(["Invalid site species"])
        # Check if all genome accession numbers come from the same taxon
        self.check_genome_accession_origin()
        # Check if all TF accession numbers come from the same taxon
        self.check_TF_accession_origin()
        return cd

class TechniquesForm(forms.Form):
    """Form to enter experimental techniques used to identify TFBS.

    - Here, curators are asked to specify ALL techniques that are used to
      identify the sites that they plan to report in the curation. The total
      number of different techniques will be used to create and populate the
      columns in the final form.
    - Curators are also asked to provide a brief description of the experimental
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
    - Curators, as before, will be prompted to specify any external DBs where
      high-throughput data might be stored (e.g. array data on GEO).
    - Curators will also be given the option to specify whether the TF is shown
      to interact with another protein/ligand that influences binding (and
      optionally add notes on that [pop-up]).
      """

    def __init__(self, *args, **kwargs):
        """Override initialization"""
        super(TechniquesForm, self).__init__(*args, **kwargs)
        help_dict = help_texts.techniques_form
        num_external_db_fields = settings.NUMBER_OF_EXTERNAL_DATABASE_FIELDS
        # Extra external-database fields
        external_db_type_choices = [(None, "None"),]
        for db in ExternalDatabase.objects.all():
            external_db_type_choices.append((db.ext_database_id,
                                             db.ext_database_name))
        # Create extra fields
        for i in xrange(num_external_db_fields):
            self.fields['external_db_type_%d' % i] = \
                forms.ChoiceField(choices=external_db_type_choices,
                                  required=False,
                                  label="External DB type [%d]" % (i+1),
                                  help_text=help_dict['external_db_type'])
            self.fields['external_db_accession_%d' % i] = \
                forms.CharField(required=False,
                                label="External DB accession number [%d]" % (i+1),
                                help_text=help_dict['external_db_accession'])

    help_dict = help_texts.techniques_form

    # Generate techniques field by getting available techniques from db
    template = get_template('experimental_technique_field.html')
    choices = [(t.technique_id,
                template.render(Context({
                    'experiment_id': t.technique_id,
                    'experiment_name': t.name,
                    'experiment_description': t.description,
                })))
               for t in ExperimentalTechnique.objects.order_by('name')]
    techniques = forms.MultipleChoiceField(
        choices = choices, label = "Techniques",
        help_text=help_dict['techniques'],
        widget=forms.CheckboxSelectMultiple())

    experimental_process = forms.CharField(widget=forms.Textarea,
                                           label="Experimental process",
                                           help_text=help_dict['experimental_process'])

    # Does TF interact with any other protein/ligand that influences binding?
    forms_complex = forms.BooleanField(required=False,
                                       label="""The manuscript reports that TF
                                       forms complex with other proteins for
                                       binding with reported sites""")

    complex_notes = forms.CharField(widget=forms.Textarea, required=False,
                                    label="Notes",
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
    help_dict = help_texts.site_entry_form
    # Type of sites to be entered
    # The curator is able to choose one of the available motifs or create a new
    # one. Sites can also be curated as either variable-motif-associated or
    # non-motif-associated.
    # To be populated dynamically
    site_type = forms.ChoiceField(choices=(), widget=forms.RadioSelect,
                                  required=True,
                                  label="Site type",
                                  help_text=mark_safe(help_dict['site_type']))

    sites = forms.CharField(required=True,
                            widget=forms.Textarea,
                            label="Sites",
                            help_text=mark_safe(help_dict['sites']))

    quantitative_data_format = forms.CharField(
        required=False,
        label="Quantitative data format",
        help_text=help_dict['quantitative_data_format'])

    # Following fields will be visible only if the submission is high-throughput
    peaks = forms.CharField(widget=forms.Textarea,
                            label="High-throughput sequences",
                            help_text=help_dict['peaks'])
    assay_conditions = forms.CharField(label="Assay conditions",
                                       help_text=help_dict['assay_conditions'],
                                       widget=forms.Textarea)
    method_notes = forms.CharField(label="Method notes",
                                   help_text=help_dict['method_notes'],
                                   widget=forms.Textarea)
    peak_techniques = forms.MultipleChoiceField(
        label="Techniques used to identify high-throughput data",
        help_text=help_dict['peak_techniques'],
        choices=(),
        widget=forms.CheckboxSelectMultiple)

    def verify_coordinates(self, coor_a, coor_b):
        try:
            x, y = int(coor_a), int(coor_b)
            if x <= 0 or y <= 0:
                return None
            return x, y
        except ValueError:
            return None

    def verify_float(self, s):
        """Given string s, check if it is in float format."""
        try:
            float(s)
            return True
        except ValueError:
            return None

    def verify_only_sites(self, sites_cd):
        """Clean function to check if all site sequences are valid."""
        try:
            if sites_cd.startswith('>'): # check if it is fasta format
                site_entry.parse_fasta(sites_cd)
            else:
                lines = [line.split() for line in sites_cd.split('\n')]
                site_entry.parse_seq(sites_cd)
        except:
            msg = """Ambiguous DNA sequence."""
            raise forms.ValidationError(msg)
        return sites_cd

    def verify_only_coordinates(self, sites_cd):
        """Clean function to check if all coordinates are valid"""
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
        """Verify input fields which have sequences and quantitative values (one
        per site)"""
        lines = [line.split() for line in sites_cd.split('\n')]
        if (all(len(line) == 2 for line in lines) and
            all(nuc in 'ACTG' for line in lines for nuc in line[0]) and
            all(self.verify_float(line[1]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("Invalid format")

        return sites_cd

    def verify_only_coordinates_and_values(self, sites_cd):
        """Verify input fields which may have coordinates and quantitative
        values (one per site)"""
        lines = [line.split() for line in sites_cd.split('\n')]
        if (all(len(line)==3 for line in lines) and
            all(self.verify_coordinates(line[0], line[1]) for line in lines) and
            all(self.verify_float(line[2]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("Invalid input format.")
        return sites_cd

    def check_quantitative_data_format(self, qval_data_format):
        """Check if the quantitative-data-format is filled."""
        if not qval_data_format:
            msg = """The reported sites have associated quantitative
            values. Please enter the quantitative data format."""
            self._errors["__all__"] = self.error_class([msg])

    def clean_sites(self):
        """Validate sites field"""
        cd = self.cleaned_data['sites'].strip().upper()

        lines = [re.split('[\t ]+', line.strip())
                 for line in re.split('[\r\n]+', cd.strip())]
        sites_cd = '\n'.join(' '.join(wds for wds in line) for line in lines)
        # by default, don't check the qval data format
        self.check_qval_data_format = False
        # check if it is fasta
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

        # else return validation error
        raise forms.ValidationError("Invalid format.")

    def clean_peaks(self):
        """Validate peaks field"""
        cd = self.cleaned_data['peaks'].strip().upper()
        lines = [re.split('[\t ]+', line.strip())
                 for line in re.split('[\r\n]+', cd.strip())]
        peaks_cd = '\n'.join(' '.join(wds for wds in line) for line in lines)
        # by default, don't check the qval data format
        self.check_qval_data_format = False
        if len(lines[0]) == 1: # sequence format
            # Same process with site verification
            return self.verify_only_sites(peaks_cd)
        elif len(lines[0]) == 2:
            # either coordinate format or sequence with qvals
            if lines[0][0].isalpha(): # seq with qvals
                self.check_qval_data_format = True
                return self.verify_only_sites_and_values(peaks_cd)
            else:
                return self.verify_only_coordinates(peaks_cd)
        elif len(lines[0]) == 3: # coordinates with qval
            self.check_qval_data_format = True
            return self.verify_only_coordinates_and_values(peaks_cd)
        else:
            raise forms.ValidationError("Invalid format")

    def clean(self):
        # If the site field contains quantitative data, make sure the format
        # field is filled.
        if 'sites' in self.cleaned_data and self.check_qval_data_format:
            self.check_quantitative_data_format(
                self.cleaned_data.get('quantitative_data_format', None))
        return self.cleaned_data

class SiteExactMatchForm(forms.Form):
    """Form to select and match reported sites to their equivalents in the
    genome. This form displays only exact matches (i.e. ones that is present in
    genome exactly). If the site is not found in the genome 'exactly', it is
    searched 'softly' and presented in the next form, SiteSoftMatchForm."""
    # all form fields are created dynamically. No static field def here
    def clean(self):
        print self.cleaned_data
        return self.cleaned_data

class SiteSoftMatchForm(forms.Form):
    """Form displaying results of 'soft' search. Match sites are not exactly
    same with the query sequence, but similar."""
    # No static def either here.
    pass

class SiteAnnotationForm(forms.Form):
    """In this form, the user is asked to fill the information regarding each
    site. In particular, the user can:
    - visualize again site information, including chromosome,
    - toggle the graphical view (off by default)
    - edit the qualitative values
    - specify which techniques were used to determine this site
    - define the site as repressed, activated or whether the effect of TF on the
      site is unknown (not-determined)

    Each site and its set of fields that can be edited are represented as a form
    (one form per site instance). The abstraction to work multiple forms in one
    page is achieved via Django FormSets
    (https://docs.djangoproject.com/en/1.6/topics/forms/formsets/).
    """
    def clean(self):
        # TODO Check if at least one experimental techniques is selected for
        # each site.
        return self.cleaned_data


class GeneRegulationForm(forms.Form):
    """This form is displayed after SiteSoftMatchForm. After the user selects
    site equivalent for each reported site, in this form, surrounding genes to
    the site are displayed. For each gene, it can be (un)checked whether site
    regualates gene (or not).

    Like the previous two forms (SiteExactMatchForm and SiteSoftMatchForm), all
    fields in this form are created dynamically, based on which genome positions
    are selected in the previous two forms as site equivalents."""
    pass


class CurationReviewForm(forms.Form):
    """Form to review all the data entered so far. The last step to submit
    curation."""
    help_dict = help_texts.curation_review_form

    choices = ((None, "None"),) + Curation.REVISION_REASONS
    revision_reasons = forms.ChoiceField(choices=choices,
                                         label="Revision required",
                                         help_text=help_dict['revision_reasons'])

    confidence = forms.BooleanField(required=False,
                                    label="I am confident of the results reported in this manuscript.",
                                    help_text=help_dict['confidence'])

    paper_complete = forms.BooleanField(required=False,
                                        label="Curation for this paper is complete.",
                                        help_text=help_dict['paper_complete'])

    notes = forms.CharField(widget=forms.Textarea,
                            required=False,
                            label="Notes",
                            help_text=help_dict['notes'])

    confirm = forms.BooleanField(required=True,
                                 label="I want to submit this curation",
                                 help_text=help_dict['confirm'])
