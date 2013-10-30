from django import forms
from collectfapp import curationform
from baseapp.templatetags import pretty_print, gene_diagram

class EditCurationForm(curationform.GenomeForm,
                       curationform.TechniquesForm,
                       curationform.CurationReviewForm):

    def __init__(self, *args, **kwargs):
        curation = kwargs.pop('curation', None)
        super(EditCurationForm, self).__init__(*args, **kwargs)
        # populate site instances
        self._populate_site_instances(curation)
        # make genome text readonly
        self.fields['genome_accession'].widget.attrs['readonly'] = True
        # hide some fields
        self.fields['site_species_same'].widget = forms.HiddenInput()
        self.fields['TF_species_same'].widget = forms.HiddenInput()

    def _populate_site_instances(self, curation):
        for csi in curation.curation_siteinstance_set.all():
            site_instance = csi.site_instance
            seq = csi.site_instance.seq
            strand = '+' if site_instance.strand == 1 else '-'
            loc = '[%d,%d]' % (site_instance.start+1, site_instance.end+1)
            label = pretty_print.site2label(csi.pk, seq+' '+strand+loc)
            help_text = gene_diagram.regulation_diagram(csi.regulation_set.all(), csi.site_instance)
            self.fields["site_instance_%d"%csi.pk] = forms.BooleanField(label=label,
                                                                        help_text=help_text,
                                                                        required=False)
            self.fields["site_instance_%d"%csi.pk].initial = True
            
    def clean(self):
        cleaned_data = super(EditCurationForm, self).clean()
        return cleaned_data

        
