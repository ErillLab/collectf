from django import forms
import models

class ExportForm(forms.Form):
    TF_instances = forms.ModelChoiceField(queryset=models.TFInstance.objects.all(),
                                label="TF",
                                initial=0)

    genomes = forms.ModelChoiceField(queryset=models.Genome.objects.all(),
                                     label="genomes",
                                     initial=0)

    is_test_export = forms.BooleanField(label="tbl export for test purposes.",
                                        help_text="If checked, site_instances will not be marked as 'submitted to NCBI'.",
                                        required=False)
