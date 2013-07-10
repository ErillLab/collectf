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
    
    TF_species_same = forms.BooleanField(required=False,
                                         label="""This is the exact same
                                         strain as reported in the manuscript for the TF.""")
    
    site_species_same = forms.BooleanField(required=False,
                                           label="""This is the exact same
                                           strain as reported in the
                                           manuscript for the sites.""") 


    TF_accession = forms.CharField(label="TF accession number",
                                   help_text=help_dict['TF_accession'])

    TF_species = forms.CharField(label="Organism of origin for reported TF",
                                 required=False,
                                 help_text=help_dict['TF_species'])

    site_species = forms.CharField(label="Organism TF binding sites are reported in",
                                   required=False,
                                   help_text=help_dict['site_species'])

    
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
            if not genome_rec:
                msg = "Cannot fetch genome record from NCBI. Check accession number."
                raise forms.ValidationError(msg)
            # create Genome object
            make_genome(genome_rec)
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
        return genome.strain.name


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
    motif_associated = forms.BooleanField(required=False,
                                          initial=True,
                                          label="Reported sites are motif associated.",
                                          help_text=help_dict['motif_associated'])
    
    is_coordinate = forms.BooleanField(initial=False,
                                       required=False,
                                       label="Binding sequence coordinates are reported in the paper.",
                                       help_text=help_dict['is_coordinate'])

    is_chip_seq_data = forms.BooleanField(initial=False,
                                          required=False,
                                          label="This paper reports Chip-Seq data",
                                          help_text=help_dict['is_chip_seq_data'])
    
    # fields about chip-seq data. To be used only if data is chip-seq
    peak_calling_method = forms.CharField(required=False,
                                          label="Peak calling method")
    
    assay_conditions = forms.CharField(required=False,
                                       label="Assay conditions",
                                       widget=forms.Textarea)
    
    method_notes = forms.CharField(required=False,
                                   label="Method notes",
                                   widget=forms.Textarea)
    
    sites = forms.CharField(required=True,
                            widget=forms.Textarea,
                            label="Sites") # fill help text with js

    def clean(self):
        cleaned_data = super(SiteReportForm, self).clean()
        if "sites" not in cleaned_data:
            return  # raise form validation error

        print cleaned_data
        c_motif_associated = cleaned_data.get("motif_associated")
        c_is_coordinate = cleaned_data.get("is_coordinate")
        c_is_chip_seq_data = cleaned_data.get("is_chip_seq_data")
        
        if c_motif_associated:
            # Reported sites are motif associated. In other words, they are true
            # binding sites that TF binds, not a arbitrarily long sequence (that have
            # the actual binding site somewhere in it) identified using Chip-Seq or
            # some other experimental method.
            sites_cd = cleaned_data.get("sites")
            sites = sitesearch.parse_site_input(sites_cd)
            if not sites:
                raise forms.ValidationError("ambiguous DNA sequence")
        else:
            # means the site is not a true site, it can be a few 100s bp and have site somewhere in it.
            if not c_is_coordinate: # reported text is a list of sequences
                sites_cd = cleaned_data.get("sites")
                sites = sitesearch.parse_site_input(sites_cd)
                if not sites:
                    raise forms.ValidationError("ambiguous DNA sequence")
            else:  # reported text is a list of coordinates for each sequences.
                sites_cd = cleaned_data.get("sites")
                while sites_cd[-1] in "\r\n\t ":
                    sites_cd = sites_cd[:-1]
                coordinates = [re.split("[\t ,]+", line) for line in re.split("[\r\n]+", sites_cd)]
                print coordinates
                if (c_is_chip_seq_data and any(len(coordinates[i]) != 3 for i in xrange(len(coordinates)))):
                    raise forms.ValidationError("Invalid coordinate input, all instances should have 3 fields (start, end, intensity-value).")
                if (not c_is_chip_seq_data and any(len(coordinates[i]) != 2 for i in xrange(len(coordinates)))):
                    raise forms.ValidationError("Invalid coordinate input, all instances should have 2 fields (start, end).")

                if any(int(coor[0]) >= int(coor[1]) for coor in coordinates):
                    raise forms.ValidationError("Invalid coordinates")

        if c_is_chip_seq_data:
            c_peak_calling_method = cleaned_data.get("peak_calling_method")
            if not c_peak_calling_method:
                msg = "Peak calling method can not be blank."
                self._errors["peak_calling_method"] = self.error_class([msg])
            c_assay_conditions = cleaned_data.get("assay_conditions")
            if not c_assay_conditions:
                msg = "Assay conditions field can not be blank."
                self._errors["assay_conditions"] = self.error_class([msg])
            c_method_notes = cleaned_data.get("method_notes")
            if not c_method_notes:
                msg = "Method notes field can not be blank."
                self._errors["method_notes"] = self.error_class([msg])
                
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
