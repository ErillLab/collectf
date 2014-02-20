from django import forms
import add_curation_form
from baseapp.templatetags import gene_diagram
from curate import site_entry

class EditCurationForm(add_curation_form.GenomeForm,
                       add_curation_form.TechniquesForm,
                       add_curation_form.CurationReviewForm):

    def __init__(self, *args, **kwargs):
        curation = kwargs.pop('curation', None)
        super(EditCurationForm, self).__init__(*args, **kwargs)
        # populate site instances
        self._populate_site_instances(curation)
        # if chip-data, add fields
        #self._populate_chip_fields(curation)
        # make genome text readonly
        self.fields['genome_accession'].widget.attrs['readonly'] = True
        # hide some fields
        self.fields['site_species_same'].widget = forms.HiddenInput()
        self.fields['TF_species_same'].widget = forms.HiddenInput()

    def _populate_site_instances(self, curation):
        for csi in curation.curation_siteinstance_set.all():
            site_instance = csi.site_instance
            # Use match class to draw the diagram
            m = site_entry.Match(csi.site_instance.genome, csi.site_instance.seq, csi.annotated_seq,
                                 csi.site_instance.start, csi.site_instance.end, csi.site_instance.strand)

            label = m.pprint()
            self.fields["site_instance_%d"%csi.pk] = forms.BooleanField(label=label,
                                                                        help_text=help_text,
                                                                        required=False)
            self.fields["site_instance_%d"%csi.pk].initial = True
            
    def clean(self):
        cleaned_data = super(EditCurationForm, self).clean()
        return cleaned_data

        
