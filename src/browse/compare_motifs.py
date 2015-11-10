"""Module for motif comparison."""

import matplotlib.pyplot as plt
import StringIO
from base64 import b64encode

from django.contrib import messages
from django.shortcuts import render
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect

from .search import get_all_TF_families
from .search import get_all_phyla
from .search import get_all_techniques
from .search import search_post_helper

from .motif_report import build_motif_reports
from .motif_report import build_ensemble_report


def motif_comparison_get(request, step):
    """Renders the form for motif comparison page."""
    form_titles = ["Step 1 of 3: Search for the first motif.",
                   "Step 2 of 3: Search for the second motif."]
    form_descriptions = [
        """Search in CollecTF is fully customizable. Just select a taxonomic
        unit (e.g. the Vibrio genus), a transcription factor family or instance
        (e.g. LexA) and the set of experimental techniques that reported sites
        should be backed by and proceed.""",
        """Search for the second motif to be compared. Just select a taxonomic
        unit (e.g. the Vibrio genus), a transcription factor family or instance
        (e.g. LexA) and the set of experimental techniques that reported sites
        should be backed by and proceed."""]
    return render(request, 'compare_motifs.html',
                  {'TF_families': get_all_TF_families(),
                   'phyla': get_all_phyla(),
                   'all_techniques': get_all_techniques(),
                   'form_title': form_titles[step-1],
                   'form_description': form_descriptions[step-1],
                   'compare': True})


def motif_comparison_search_first_motif(request):
    """View for step 1 of the motif comparison: search the first motif."""
    if request.POST:
        return motif_comparison_post(request, 1)
    return motif_comparison_get(request, 1)


def motif_comparison_search_second_motif(request):
    """View for step 2 of the motif comparison: search the second motif"""
    if request.POST:
        return motif_comparison_post(request, 2)
    if 'motif_comparison_step_1' not in request.session:
        return motif_comparison_get(request, 1)
    return motif_comparison_get(request, 2)


def motif_comparison_post(request, step):
    """Processes motif comparison forms and compares motifs."""
    try:
        curation_site_instances = search_post_helper(request)
    except ValidationError, e:
        messages.add_message(request, messages.ERROR, e.message)
        return motif_comparison_get(request, step)

    if not curation_site_instances.filter(site_type='motif_associated'):
        message = """No results found for selected TF, species and experimental
        techniques."""
        messages.add_message(request, messages.ERROR, message)
        return motif_comparison_get(request, step)

    # Store search results
    request.session['motif_comparison_step_%d' % step] = curation_site_instances
    request.session.modified = True

    if step == 1:
        return HttpResponseRedirect(reverse(
            motif_comparison_search_second_motif))
    else:
        # Step 2
        # Render results
        first_sites = request.session['motif_comparison_step_1']
        second_sites = request.session['motif_comparison_step_2']

        return render(
            request, 'compare_motifs_results.html',
            {'form_title': "Step 3 of 3: Motif comparison results",
             'first_motif_reports': build_motif_reports(first_sites),
             'first_ensemble_report': build_ensemble_report(first_sites),
             'second_motif_reports': build_motif_reports(second_sites),
             'second_ensemble_report': build_ensemble_report(second_sites)})


def motif_similarity_measure(request):
    """AJAX handler for motif similarity measurements.

    Request has the similarity function and two motifs (space separated
    strings). Returns HTML which will be embedded into the main comparison
    page.
    """

    unaligned_sites_first = request.POST['first_unaligned'].strip().split()
    unaligned_sites_second = request.POST['second_unaligned'].strip().split()
    sites_first = request.POST['first_aligned'].strip().split()
    sites_second = request.POST['second_aligned'].strip().split()

    print unaligned_sites_first

    if request.POST['similarity_function'] == 'site_based':
        boxplot, histogram = levenshtein_measure(sites_first, sites_second)
        return render(request, 'motif_similarity_levenshtein.html',
                      {'boxplot': boxplot,
                       'histogram': hist,
                       'first_unaligned': unaligned_sites_first,
                       'second_unaligned': unaligned_sites_second,
                       'first_aligned': sites_first,
                       'second_aligned': sites_second})

    else:
        # other similarity metrics
        fun_str = request.POST['similarity_function']
        fun2call_dict = {
            'euclidean': euclidean_distance,
            'pearson': pearson_correlation_coefficient,
            'kullback_leibler': kullback_leibler_divergence,
            'average_log_likelihood': average_log_likelihood_ratio}
        fun2call = fun2call_dict[fun_str]
        sites_a, sites_b = motif_alignment(sites_a, sites_b)
        fig, p_val = motif_sim_test(request, sites_a, sites_b, fun2call)
        return render(request, 'motif_similarity_%s.html' % fun_str,
                      {'plot': fig,
                       'sc': '%.3lf' % fun2call(sites_a, sites_b),
                       'p_value': '%.3lf' % p_val,
                       'unaligned_sites_a': unaligned_sites_a,
                       'unaligned_sites_b': unaligned_sites_b,
                       'aligned_sites_a': sites_a,
                       'aligned_sites_b': sites_b})


def levenshtein_distance(motif_a, motif_b):
    """
    Given two sets of binding sites (motif_a and motif_b), measure Levenshtein
    distance between all pairs of sequences from two motifs, to measure the distance
    between two motifs
    """
    dists = [float(levenshtein(seqa, seqb)) / (len(seqa) * len(seqb))
            for seqa in motif_a for seqb in motif_b]
    return dists

def levenshtein_measure(motif_a, motif_b):
    """Levenshtein distance between two motifs."""
    a_vs_a = levenshtein_distance(motif_a, motif_a)
    a_vs_b = levenshtein_distance(motif_a, motif_b)
    b_vs_b = levenshtein_distance(motif_b, motif_b)

    bp = plt.boxplot([a_vs_a, a_vs_b, b_vs_b])
    plt.xticks(range(1, 4),
               [r'$M_a$ vs $M_a$', r'$M_a$ vs $M_b$', r'$M_b$ vs $M_b$'])
    plt.setp(bp['boxes'], color='#0088CC')
    plt.setp(bp['whiskers'], color='#0088CC')
    plt.setp(bp['fliers'], color='#0088CC', marker='+')
    plt.ylabel('Levenshtein distance')
    boxplot = fig2img(plt.gcf())
    plt.hist([a_vs_a, a_vs_b, b_vs_b], bins=15,
             label=[r'$M_a$ vs $M_a$', r'$M_a$ vs $M_b$', r'$M_b$ vs $M_b$'])
    plt.ylabel('Frequency of site-pairs')
    plt.xlabel('Levenshtein distance')
    plt.legend()
    hist = fig2img(plt.gcf())
    return boxplot, hist


def fig2img(fig):
    """Returns data URI, given the matplotlib plot."""
    imgdata = StringIO.StringIO()
    fig.savefig(imgdata, format='png')
    plt.clf()
    imgdata.seek(0) # rewind the data
    return "data:image/png;base64,%s" % b64encode(imgdata.buf)
