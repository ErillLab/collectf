"""Module for scanning a given genome with a motif."""

from django.shortcuts import render
from django.views.decorators.http import require_POST

from core import bioutils

from .forms.pssm_search import BindingSiteSearchForm


@require_POST
def pssm_search_from_report_page(request):
    """Gets the list of sites and renders the PSSM search form with populated
    sites field."""
    sites = request.POST['sites'].split()
    form = BindingSiteSearchForm({'sites': '\n'.join(sites)})

    return render(request, 'pssm_search.html', {'form': form})


def pssm_search(request):
    """Scans the genome with the given motif."""
    if request.method == 'POST':
        form = BindingSiteSearchForm(request.POST)
        if form.is_valid():
            genome = form.cleaned_data['genome']
            motif = bioutils.build_motif(form.cleaned_data['sites'])
            # Find a threshold
            dist = motif.pssm.distribution(precision=10**4)
            threshold = dist.threshold_patser()
            hits = bioutils.pssm_search(
                motif.pssm, genome.seq, threshold=threshold)
            matches = []
            for pos, strand, score in hits:
                seq = genome.seq[pos:pos+motif.length]
                if strand == -1:
                    seq = bioutils.reverse_complement(seq)
                matches.append({
                    'start': pos,
                    'end': pos + motif.length,
                    'strand': strand,
                    'seq': seq,
                    'score': score})
            return render(
                request,
                'pssm_search_results.html',
                {'matches': matches,
                 'weblogo': bioutils.weblogo_uri(form.cleaned_data['sites']),
                 'threshold': threshold})
    else:
        form = BindingSiteSearchForm()

    return render(request, 'pssm_search.html', {'form': form})
