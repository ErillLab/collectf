# Forms for curation wizard
from django import forms
from models import *
import bioutils
from makeobj import *


class PublicationForm(forms.Form):
    """Publication selection form"""
    pub = forms.ChoiceField(widget=forms.RadioSelect())

class GenomeForm(forms.Form):
    """Form for submission of genome and TF accession numbers and others"""
    genome_accession = forms.CharField(label="genome accession")
    TF_accession = forms.CharField(label="TF accession")
    TF_function = forms.ChoiceField(Curation.TF_FUNCTION, label="TF function")
    TF_type = forms.ChoiceField(Curation.TF_TYPE, label="TF type")
    TF = forms.ChoiceField(choices=(), label="TF")
    
    # form fields below are redundant.  When the form is displayed, the curator
    # can either input TF and site species manually, or he can select the option
    # that they are same species with genome. If they are same species with
    # genome, he doesn't need to input whole species name
    TF_species_same = forms.BooleanField(label="TF-species is same with genome",
                                         required=False)
    TF_species = forms.CharField(label="TF species", required=False)
    # neither is required, but one of them should be filled
    site_species_same = forms.BooleanField(label="site-species is same with genome",
                                           required=False)
    site_species = forms.CharField(label="site species", required=False)

    def clean_genome_accession(self):
        genome_accession = self.cleaned_data['genome_accession']
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
        TF_accession = self.cleaned_data['TF_accession']
        try:
            TF_instance = TFInstance.objects.get(protein_accession=TF_accession)
        except TFInstance.DoesNotExist:
            TF_rec = bioutils.get_TF(TF_accession)
            if not TF_rec:
                msg = "Cannot fetch protein record from NCBI. Check accession number."
                raise forms.ValidationError(msg)
            # create TFInstance object
        return TF_accession
