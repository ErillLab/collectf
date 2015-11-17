from django import forms
from django.template import Context
from django.template.loader import get_template

from . import help_text
from core.models import ExperimentalTechnique, ExternalDatabase

class TechniquesForm(forms.Form):
    """Form to enter experimental techniques used to identify TFBS.
    
    Curators are asked to specify ALL techniques that are used to identify the
    sites that they plan to report in the curation. 
    
    Curators are also asked to provide a brief description of the experimental
    setup for the sites reported. The description should comprise the setup
    used for all the sites that are going to be reported, even though they
    might used different techniques. For instance: "The binding motif for XX
    was identified through phylogenetic footprinting and directed-mutagenesis +
    EMSA on the promoter of gene YY. Researchers then performed a computer
    search of sites with one mismatch in the genome. Of those identified, they
    verified X through EMSA. Regulatory activity for four of the sites was
    assessed with beta-gal assays.". In brief, the curator is expected to
    provide a concise and logical summary of the experimental process leading
    to the identification of all reported sites and the determination (if any)
    of their regulatory activity.  Curators, as before, will be prompted to
    specify any external DBs where high-throughput data might be stored
    (e.g. array data on GEO).
    
    Curators will also be given the option to specify whether the TF is shown
    to interact with another protein/ligand that influences binding (and
    optionally add notes on that [pop-up]).
    """

    NUM_EXTERNAL_DB_FIELDS = 5

    def __init__(self, *args, **kwargs):
        """Overrides initialization."""        
        super(TechniquesForm, self).__init__(*args, **kwargs)

        # Extra external-database fields
        external_db_type_choices = [(None, 'None'),]
        for db in ExternalDatabase.objects.all():
            external_db_type_choices.append(
                (db.ext_database_id, db.ext_database_name))
        # Create extra fields
        for i in range(self.NUM_EXTERNAL_DB_FIELDS):
            self.fields['external_db_type_%d' % i] = forms.ChoiceField(
                choices=external_db_type_choices,
                required=False,
                label="External DB type (%d)" % i,
                help_text=TechniquesForm.help_text_dict['external_db_type'])
            self.fields['external_db_accession_%d' % i] = forms.CharField(
                required=False,
                label="External DB accession number (%d)" % i,
                help_text=TechniquesForm.help_text_dict['external_db_accession'])

    help_text_dict = help_text.techniques_form

    # Generate techniques field by getting available techniques from db
    template = get_template('experimental_technique_field.html')
    choices = [(technique.technique_id,
                template.render(Context({'technique': technique})))
               for technique in ExperimentalTechnique.objects.order_by('name')]
    
    techniques = forms.MultipleChoiceField(
        choices=choices,
        label="Techniques",
        help_text=help_text_dict['techniques'],
        widget=forms.CheckboxSelectMultiple())

    experimental_process = forms.CharField(
        widget=forms.Textarea,
        label="Experimental process",
        help_text=help_text_dict['experimental_process'])

    # Does TF interact with any other protein/ligand that influences binding?
    forms_complex = forms.BooleanField(
        required=False,
        label="""
        The manuscript reports that TF forms complex with other proteins for
        binding with reported sites.""")

    complex_notes = forms.CharField(
        widget=forms.Textarea, required=False, label="Notes",
        help_text=help_text_dict['complex_notes'])

    # External database links
    has_external_db = forms.BooleanField(
        required=False,
        label="""The manuscript reports high-throughput data from an external
        database. (You can report up to 5 external resources.)""")

