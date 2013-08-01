# Forms for curation wizard
from django import forms
from models import *
import bioutils
from makeobj import *
import sitesearch
from django.utils.safestring import mark_safe
import re
import helptext

class PublicationForm(forms.Form):
    """Publication selection form"""
    help_dict = helptext.publication_form
    pub = forms.ChoiceField(widget=forms.RadioSelect(),
                            label="Publications",
                            help_text=help_dict['pub'])
    
    no_data = forms.BooleanField(label="This paper contains no data.",
                                 required=False,
                                 help_text=help_dict['no_data'])

class GenomeForm(forms.Form):
    """Form for submission of genome and TF accession numbers and others"""
    help_dict = helptext.genome_form
    TF = forms.ModelChoiceField(queryset=TF.objects.order_by('name'),
                                label="TF",
                                help_text=help_dict['TF'])
    
    TF_type = forms.ChoiceField(Curation.TF_TYPE,
                                label="TF structure",
                                help_text=help_dict['TF_type'])
    
    TF_function = forms.ChoiceField(Curation.TF_FUNCTION,
                                    label="TF function",
                                    help_text=help_dict['TF_function'])

    genome_accession = forms.CharField(label="Genome NCBI accession number",
                                       help_text=help_dict['genome_accession'])

    site_species_same = forms.BooleanField(required=False,
                                           label="""This is the exact same
                                           strain as reported in the
                                           manuscript for the sites.""")   

    TF_accession = forms.CharField(label="TF accession number",
                                   help_text=help_dict['TF_accession'])

    TF_species_same = forms.BooleanField(required=False,
                                         label="""This is the exact same
                                         strain as reported in the manuscript for the TF.""")
    
    site_species = forms.CharField(label="Organism TF binding sites are reported in",
                                   required=False,
                                   help_text=help_dict['site_species'])   

    TF_species = forms.CharField(label="Organism of origin for reported TF",
                                 required=False,
                                 help_text=help_dict['TF_species'])


    
    def clean_genome_accession(self):
        """Validate entered genome accession number."""
        genome_accession = self.cleaned_data['genome_accession'].strip()
        
        if '.' not in genome_accession:
            msg = 'Please enter RefSeq accession number with version number'
            raise forms.ValidationError(msg)
        
        try: # to retrieve genome from database (if exists)
            g = Genome.objects.get(genome_accession=genome_accession)
        except Genome.DoesNotExist: # try to retrieve from NCBI database
            genome_rec = bioutils.get_genome(genome_accession)
            strain_taxon = bioutils.get_org_taxon(genome_rec)
            if not genome_rec:
                msg = "Cannot fetch genome record from NCBI. Check accession number."
                raise forms.ValidationError(msg)
            if not strain_taxon:
                msg = "Cannot fetch strain taxonomy information."
                raise forms.ValidationError(msg)
            # create Genome object
            make_genome(genome_rec, strain_taxon)
            # create all Genome objects for each gene on that genome
            make_all_genes(genome_rec)
        return genome_accession

    def clean_TF_accession(self):
        TF_accession = self.cleaned_data['TF_accession'].strip()
        # remove version
        TF_accession = TF_accession.split('.')[0]
        try:
            TF_instance = TFInstance.objects.get(protein_accession=TF_accession)
        except TFInstance.DoesNotExist:
            TF_rec = bioutils.get_TF(TF_accession)
            if not TF_rec:
                msg = "Cannot fetch protein record from NCBI. Check accession number."
                raise forms.ValidationError(msg)
            # create TFInstance object
            make_TF_instance(TF_rec)
        return TF_accession

    def clean(self):
        # check if either site-species or site-species-same is filled
        # check if either TF-species or TF-species-same is filled
        cd = self.cleaned_data
        if not (cd['TF_species'] or cd['TF_species_same']):
            self._errors['TF_species'] = self.error_class([u"Invalid TF species"])
        if not (cd['site_species'] or cd['site_species_same']):
            self._errors['site_species'] = self.error_class([u"Invalid site species"])
        return cd

    def clean_TF_species(self):
        # clean TF_species field. If TF_species_same checkbox is selected, assign
        # species data to cleaned data. In that case, it is important that
        # clean_genome is called BEFORE clean_TF_species, because genome is
        # needed to extract species information.
        cd = self.cleaned_data
        if cd['TF_species_same']:
            tf_species = self.clean_species("TF_species")
        else:
            tf_species = cd['TF_species']
        return tf_species

    def clean_site_species(self):
        # same with clean_TF_species function, but for site_species field
        cd = self.cleaned_data
        if cd['site_species_same']:
            site_species = self.clean_species("site_species")
        else:
            site_species = cd['site_species']
        return site_species

    def clean_species(self, field):
        # Helper function for clean_TF_species and clean_site_species
        cd = self.cleaned_data
        if 'genome_accession' not in self.cleaned_data:
            return
        genome_accession = self.cleaned_data['genome_accession']
        genome = Genome.objects.get(genome_accession=genome_accession)
        if not genome:
            msg = u"Invalid genome accession number"
            self._errors[field] = self.error_class([msg])
        return genome.organism


class TechniquesForm(forms.Form):
    """Form to enter experimental techniques used to identify TFBS"""
    help_dict = helptext.techniques_form
    contains_promoter_data = forms.BooleanField(required=False,
                                                label="The manuscript contains promoter information",
                                                help_text=help_dict['contains_promoter_data'])

    contains_expression_data = forms.BooleanField(required=False,
                                                  label="The manuscript contains expression data",
                                                  help_text=help_dict['contains_expression_data'])
    # generate techniques field
    # get available techniques from db
    choices = []
    # used twitter-bootstrap tooltip for exp technique description
    description_markup = u'<span data-toggle="popover" title="%s" data-content="%s">%s</span>'
    for t in ExperimentalTechnique.objects.order_by('name'):
        choices.append((t.technique_id,
                        mark_safe(description_markup % (t.name, t.description, t.name))))
    techniques = forms.MultipleChoiceField(choices = choices,
                                           label = "Techniques",
                                           help_text=help_dict['techniques'],
                                           widget = forms.CheckboxSelectMultiple())

    experimental_process = forms.CharField(widget=forms.Textarea,
                                           required=False,
                                           label="Experimental process",
                                           help_text=help_dict['experimental_process'])

    external_db_type_choices = [(None, "None"),]
    for db in ExternalDatabase.objects.all():
        external_db_type_choices.append((db.ext_database_id, db.ext_database_name))
    external_db_type = forms.ChoiceField(choices=external_db_type_choices,
                                         required=False,
                                         label="External DB type",
                                         help_text=help_dict['external_db_type'])
    
    external_db_accession = forms.CharField(required=False,
                                           label="External DB accession number",
                                           help_text=help_dict['external_db_accession'])
    
    forms_complex = forms.BooleanField(required=False,
                                       label="""The manuscript reports that TF forms complex
                                       with other proteins for binding with reported sites""")
    
    complex_notes = forms.CharField(widget=forms.Textarea, required=False,
                                    label="Notes",
                                    help_text=help_dict['complex_notes'])


class SiteReportForm(forms.Form):
    """Form to input the list of sites reported in the paper"""
    help_dict = helptext.site_report_form
    description_markup = u'<span data-toggle="tooltip" title="%s">%s</span>'

    # for these fields, show help text as popover
    is_motif_associated = forms.BooleanField(required=False,
                                             initial=True,
                                             label=description_markup%(help_dict['is_motif_associated'],
                                                                       "Reported sites are motif associated."))
                                                                       
                                             #help_text=help_dict['is_motif_associated'])

    is_chip_data = forms.BooleanField(initial=False,
                                      required=False,
                                      label=description_markup%(help_dict['is_chip_data'],
                                                                "This paper reports ChIP data"))
                                      #label="This paper reports ChIP data",
                                      #help_text=help_dict['is_chip_data'])

    has_quantitative_data = forms.BooleanField(initial=False,
                                               required=False,
                                               label=description_markup%(help_dict['has_quantitative_data'],
                                                                         "Sites with quantitative data"))
                                               #label="Sites with quantitative data",
                                               #help_text=help_dict['has_quantitative_data'])
    
    is_coordinate = forms.BooleanField(initial=False,
                                       required=False,
                                       label=description_markup%(help_dict['is_coordinate'],
                                                                 "Coordinate entry mode"))
                                       #label="Coordinate entry mode",
                                       #help_text=help_dict['is_coordinate'])

    sites = forms.CharField(required=True,
                            widget=forms.Textarea,
                            label="Sites",
                            help_text=help_dict['sites'])

    # ChIP Fields
    quantitative_data_format = forms.CharField(required=False,
                                          label=description_markup%(help_dict['quantitative_data_format'],
                                                                    "Quantitative data format"))
                                          #label="Quantitative data format")

    assay_conditions = forms.CharField(required=False,
                                       label=description_markup%(help_dict['assay_conditions'],
                                                                 "Assay conditions"),
                                       #label="Assay conditions",
                                       widget=forms.Textarea)
    
    chip_method_notes = forms.CharField(required=False,
                                        label=description_markup%(help_dict['chip_method_notes'],
                                                                  "ChIP method notes"),
                                        #label="ChIP Method notes",
                                        widget=forms.Textarea)

    # Extra ChIP-field for associating quantitative values from non-motif reported
    # ChIP sites to motif-associated ones. If user selects both ChIP-data and
    # motif-associated-site and quantitative-data options, it means that, user wants
    # to enter both sites (without peak intensities) and longer regions (with peak
    # intensities). The idea here is to associate peak_intensity values from longer
    # regions to the motif_associated ones.
    chip_data_extra_field = forms.CharField(required=False,
                                            label=description_markup%(help_dict['chip_data_extra_field'],
                                                                      "Supporting ChIP quantitative data"),
                                            #label="Supporting ChIP quantitative data.",
                                            #help_text=help_dict['chip_data_extra_field'],
                                            widget=forms.Textarea)

    # End of ChIP fields



    # Based on variables about site data, namely
    # - is_motif_associated
    # - is_chip_data
    # - has_quantitative_data
    # - is_coordinate
    # the way that form is handled properly (at least trying).
    # To avoid long and ugly if/else clauses, all cases are handled in separate
    # functions which are called from clean(self) function.

    # Before introducing those clean_helper functions, here are some verify functions
    def verify_coordinates(self, coor_a, coor_b):
        try:
            x,y = int(coor_a), int(coor_b)
            return x,y
        except ValueError:
            return None

    def verify_float(self, s):
        # given string s, check if it is in float format
        try:
            x = float(s)
            return x
        except ValueError:
            return None

    # verify_site_sequences -- in sitesearch.py (handles FASTA format too.)


    # -- helper helper functions
    def verify_only_sites(self, cleaned_data):
        # Clean function to check if all site sequences are valid
        sites_cd = cleaned_data.get('sites').strip()
        sites = sitesearch.parse_site_input(sites_cd)
        if not sites:
            msg = "Ambiguous DNA sequence"
            raise forms.ValidationError(msg)
        return cleaned_data

    def verify_only_coordinates(self, cleaned_data):
        # Clean function to check if all coordinates are valid
        sites_cd = cleaned_data.get('sites').strip()
        coordinates = [re.split('[\t ]+', line) for line in re.split('[\r\n]+', sites_cd)]
        msg = "Invalid coordinate format."
        for instance in coordinates:
            if len(instance) != 2:
                raise forms.ValidationError(msg)
            if not self.verify_coordinates(instance[0], instance[1]):
                raise forms.ValidationError(msg)
        return cleaned_data

    def verify_only_sites_and_values(self, cleaned_data):
        # Verify input fields which may have
        # - either only sequences, or
        # - sequences and quantitative values (one per site)
        sites_cd = cleaned_data.get('sites').strip()
        lines = [re.split('[\t ]+', line) for line in re.split('[\r\n]+', sites_cd)]
        if (all(len(line) == 2 for line in lines) and
            all(nuc in 'ACTGactg' for line in lines for nuc in line[0]) and
            all(self.verify_float(line[1]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("invalid format")
        
        return cleaned_data

    def verify_only_coordinates_and_values(self, cleaned_data):
        # Verify input fields which may have
        # - either only coordinates, or
        # - coordinates and quantitative values (one per site)
        sites_cd = cleaned_data.get('sites').strip()
        lines = [re.split('[\t ]+', line) for line in re.split('[\r\n]+', sites_cd)]
        if (all(len(line)==3 for line in lines) and
            all(self.verify_coordinates(line[0], line[1]) for line in lines) and
            all(self.verify_float(line[2]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("Invalid input format.")
        return cleaned_data

    def verify_chip_data_extra_field(self, cleaned_data):
        # verify chip_data extra field
        # This field is used to associate motif-associated-sites with quantitative values, 
        chip_data_extra_field = cleaned_data.get('chip_data_extra_field').strip()
        lines = [re.split('[\t ]+', line) for line in re.split('[\r\n]+', chip_data_extra_field)]
        if (all(len(line)==3 for line in lines) and
            all(self.verify_coordinates(line[0], line[1]) for line in lines) and
            all(self.verify_float(line[2]) for line in lines)):
            pass
        else:
            raise forms.ValidationError("Invalid ChIP data extra field format.")


    # real clean_helper functions (16 of them -- wow)
    def clean_helper_0(self, cleaned_data):
        # NOT is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        # Just site sequences are expected. Nothing to handle specially.
        return self.verify_only_sites(cleaned_data)
    
    def clean_helper_1(self, cleaned_data):
        # NOT is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # is_coordinate
        # The site input field contains only coordinates (start,end)
        return self.verify_only_coordinates(cleaned_data)

    def clean_helper_2(self, cleaned_data):
        # NOT is_motif_associated
        # NOT is_chip_data
        # has_quantitative_data
        # NOT is_coordinate
        # In this input type, each site is reported with sequence and quantitative data.
        # after soft_site_search (in that form, for each matched site instance, quantitative values are asked
        # to be entered.
        # Check if all fields have either only sequence or, sequence and quantitative values
        return self.verify_only_sites_and_values(cleaned_data)

    def clean_helper_3(self, cleaned_data):
        # NOT is_motif_associated
        # NOT is_chip_data
        # has_quantitative_data
        # is_coordinate
        # Coordinates may be entered in this step or later. Therefore, all lines either
        # - must have only coordinates, or
        # - must have coordinates + quantitative values
        return self.verify_only_coordinates_and_values(cleaned_data)

    def clean_helper_4(self, cleaned_data):
        # NOT is_motif_associated
        # is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        # only sequence
        return self.verify_only_sites(cleaned_data)

    def clean_helper_5(self, cleaned_data):
        # NOT is_motif_associated
        # is_chip_data
        # NOT has_quantitative_data
        # is_coordinate
        # only coordinates
        return self.verify_only_coordinates(cleaned_data)

    def clean_helper_6(self, cleaned_data):
        # NOT is_motif_associated
        # is_chip_data
        # has_quantitative_data
        # NOT is_coordinate
        # If it is ChIP but not reported with coordinates, we don't handle this case yet.
        # For now, ChIP data must be with coordinates (unless it is is_motif_associated)
        return self.verify_only_sites_and_values(cleaned_data)

    def clean_helper_7(self, cleaned_data):
        # NOT is_motif_associated
        # is_chip_data
        # has_quantitative_data
        # is_coordinate
        # Both coordinates and quantitative data. However, quantitative data may be entered later.
        return self.verify_only_coordinates_and_values(cleaned_data)

    def clean_helper_8(self, cleaned_data):
        # is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        # Just site sequences are expected. Nothing to handle specially.
        return self.verify_only_sites(cleaned_data)

    def clean_helper_9(self, cleaned_data):
        # is__motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # is_coordinate
        # ONLY coordinates are reported, nothing more.
        return self.verify_only_coordinates(cleaned_data)

    def clean_helper_10(self, cleaned_data):
        # is_motif_associated
        # NOT is_chip_data
        # has_quantitative_data
        # NOT is_coordinate
        # If this is the case, nothing special needs to be handled.
        return self.verify_only_sites_and_values(cleaned_data)

    def clean_helper_11(self, cleaned_data):
        # is_motif_associated
        # NOT is_chip_data
        # has_quantitative_data
        # is_coordinate
        return self.verify_only_coordinates_and_values(cleaned_data)

    def clean_helper_12(self, cleaned_data):
        # is_motif_associated
        # is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        return self.verify_only_sites(cleaned_data)

    def clean_helper_13(self, cleaned_data):
        # is_motif_associated
        # is_chip_data
        # NOT has_quantitative_data
        # is_coordinate
        # The site input field contains only coordinates
        return self.verify_only_coordinates(cleaned_data)

    def clean_helper_14(self, cleaned_data):
        # is_motif_associated
        # is_chip_data
        # has_quantitative_data
        # NOT is_coordinate
        self.verify_only_sites(cleaned_data)
        # extra check on ChIP extra field.
        # Always coordinates and quantitative data are expected here.
        self.verify_chip_data_extra_field(cleaned_data)
        return cleaned_data

    def clean_helper_15(self, cleaned_data):
        # is_motif_associated
        # is_chip_data
        # has_quantitative_data
        # is_coordinate
        self.verify_only_coordinates(cleaned_data)
        # extra check on ChIP extra field.
        # Always coordinates and quantitative data are expected here.
        self.verify_chip_data_extra_field(cleaned_data)
        return cleaned_data

    def clean(self):
        # Form fields that depend on each other are handled in this function.
        cleaned_data = super(SiteReportForm, self).clean()
        if "sites" not in cleaned_data:
            return  # raise form validation error
        
        cleaned_data['sites'].strip()  # remove trailing spaces
        if not cleaned_data['sites']: return # raise form validation error

        print cleaned_data

        is_motif_associated = cleaned_data.get("is_motif_associated")
        is_chip_data = cleaned_data.get("is_chip_data")
        has_quantitative_data = cleaned_data.get('has_quantitative_data')
        is_coordinate = cleaned_data.get("is_coordinate")

        # Input validation
        # Don't let user proceeed if it is not exact genome and coordinate data is selected.
        # TODO

        func = (('1' if is_motif_associated else '0') +
                ('1' if is_chip_data else '0') +
                ('1' if has_quantitative_data else '0') +
                ('1' if is_coordinate else '0'))
        func_call_str = 'self.clean_helper_%d(cleaned_data)' % int(func,2)
        print 'form_validation:', func_call_str
        eval(func_call_str)

        # extra validation for quantitative data
        if has_quantitative_data:
            if not cleaned_data.get('quantitative_data_format'):
                msg = "This field can not be blank."
                self._errors["quantitative_data_format"] = self.error_class([msg])

        # extra validations for ChIP data
        if is_chip_data:
            if not cleaned_data.get('assay_conditions'):
                msg = "Assay conditions field can not be blank."
                self._errors["assay_conditions"] = self.error_class([msg])
            if not cleaned_data.get('chip_method_notes'):
                msg = "Method notes field can not be blank."
                self._errors["chip_method_notes"] = self.error_class([msg])

                
        return self.cleaned_data

            
class SiteExactMatchForm(forms.Form):
    """Form to select and match reported sites to their equivalents in the
    genome. This form displays only exact matches (i.e. ones that is present in
    genome exactly). If the site is not found in the genome 'exactly', it is
    searched 'softly' and presented in the next form, SiteSoftMatchForm."""
    # all form fields are created dynamically. No static field def here
    pass

class SiteSoftMatchForm(forms.Form):
    """Form displaying results of 'soft' search. Match sites are not exactly
    same with the query sequence, but similar."""
    # No static def either here.
    pass

class SiteRegulationForm(forms.Form):
    """This form is displayed after SiteSoftMatchForm. After the user selects
    site equivalent for each reported site, in this form, surrounding genes to
    the site are displayed. For each gene, it can be (un)checked whether site
    regualates gene (or not).

    Like the previous two forms, all fields in this form are created
    dynamically, based on which genome positions are selected in the previous
    two forms as site equivalents."""
    pass

class CurationReviewForm(forms.Form):
    """Form to see all data entered so far. The last step to submit curation."""
    help_dict = helptext.curation_review_form
    
    choices = ((None, "None"),) + Curation.REVISION_REASONS
    revision_reasons = forms.ChoiceField(choices=choices,
                                         label="Revision required",
                                         help_text=help_dict['revision_reasons'])
    
    confidence = forms.BooleanField(required=False,
                                    label="I am confident of the results reported in this manuscript.",
                                    help_text=help_dict['confidence'])

    NCBI_submission_ready = forms.BooleanField(required=False,
                                               label="Curation is ready to submit to NCBI.",
                                               help_text=help_dict['NCBI_submission_ready'])
                                               
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


# One extra form goes here: Site quantitative value form
class SiteQuantitativeDataForm(forms.Form):
    """In this form, all matched sites are presented with prefilled quantitative
    values (if any). The purpose is to let user review quantitative data associated
    with each site, before writing anything to database. This form is displayed only
    if the curation is marked as 'sites with quantitative data'.
    """
    pass
    
