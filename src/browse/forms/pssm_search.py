from django import forms
from django.utils.safestring import mark_safe

from core import entrez_utils


class BindingSiteSearchForm(forms.Form):
    sites = forms.CharField(
        label="Sites",
        widget=forms.Textarea(attrs={'class': 'sequence'}))

    genome = forms.CharField(
        label='NCBI RefSeq accession number',
        max_length=20,
        help_text=mark_safe("""
        Enter the RefSeq accession number of the sequence to be scanned for
        binding sites, ex: <code>NC_000913</code>. <a
        href='http://www.ncbi.nlm.nih.gov/genome/browse/' target='blank_'> Find
        a genome to search.</a>"""))

    def clean_sites(self):
        input_field = self.cleaned_data['sites'].upper().strip()
        motif_sites = [site.strip() for site in input_field.split('\n')]
        # Check each line contains a valid DNA sequence.
        for site in motif_sites:
            if any(b not in 'ACGT' for b in site):
                raise forms.ValidationError(
                    "Input motif must be DNA sequences")

        # Check lengths of sequences
        if any(len(site) != len(motif_sites[0]) for site in motif_sites):
            raise forms.ValidationError(
                "Input sequences should have same length.")
        return motif_sites

    def clean_genome(self):
        try:
            genome = entrez_utils.get_genome(
                self.cleaned_data['genome'].strip())
        except entrez_utils.EntrezException:
            raise forms.ValidationError("Invalid RefSeq accession number.")
        return genome
