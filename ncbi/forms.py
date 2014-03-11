from django import forms
from base import models

class ExportForm(forms.Form):
    genome_accession = forms.CharField(label="Genome accession", required=True)

    is_test_export = forms.BooleanField(label="tbl export for test purposes.",
                                        help_text="If checked, site_instances will not be marked as 'submitted to NCBI'.",
                                        initial=True,
                                        required=False)

    def clean_genome_accession(self):
        """Check genome accession field"""
        genome_accession = self.cleaned_data['genome_accession'].strip()
        try:
            g = models.Genome.objects.get(genome_accession=genome_accession)
        except models.Genome.DoesNotExist:
            msg = "The genome with accession number %s is not in the database." % genome_accession
            raise forms.ValidationError(msg)
        return genome_accession

    def clean(self):
        return self.cleaned_data
