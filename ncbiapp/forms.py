from django import forms
import models

class ExportForm(forms.Form):
    genomes = forms.ModelChoiceField(queryset=models.Genome.objects.all(),
                                     label="Genome",
                                     initial=0)

    is_test_export = forms.BooleanField(label="tbl export for test purposes.",
                                        help_text="If checked, site_instances will not be marked as 'submitted to NCBI'.",
                                        initial=True,
                                        required=False)
