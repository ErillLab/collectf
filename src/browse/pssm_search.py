"""Module for scanning a given genome with a motif."""

from django.shortcuts import render

from core import bioutils

from .forms.pssm_search import BindingSiteSearchForm
        
def pssm_search(request):

    """Scans the genome with the given motif."""
    if request.method == 'POST':
        form = BindingSiteSearchForm(request.POST)
        if form.is_valid():            
            genome = form.cleaned_data['genome']
            motif = bioutils.build_motif(form.cleaned_data['sites'])
            # Scan the genome for binding sites.
            for pos, score in bioutils.pssm_search(
                    motif.pssm, str(genome.seq), threshold=5.0):
                print pos, score
            
    else:
        form = BindingSiteSearchForm()

    return render(request, 'pssm_search.html', {'form': form})

