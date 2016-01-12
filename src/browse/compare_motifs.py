"""Module for motif comparison."""

from base64 import b64encode
import StringIO
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random

from django.contrib import messages
from django.shortcuts import render
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.shortcuts import redirect

from core import bioutils
from .search_motifs import get_all_TF_families
from .search_motifs import get_all_phyla
from .search_motifs import get_all_techniques
from .search_motifs import search_post_helper
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
    request.session['motif_comparison_step_%d' % step] = (
        curation_site_instances)
    request.session.modified = True

    if step == 1:
        return redirect(motif_comparison_search_second_motif)
    else:
        # Step 2
        # Render results
        first_sites = request.session['motif_comparison_step_1']
        second_sites = request.session['motif_comparison_step_2']

        first_motif_reports = build_motif_reports(first_sites)
        second_motif_reports = build_motif_reports(second_sites)
        first_ensemble_report = build_ensemble_report(first_sites)
        second_ensemble_report = build_ensemble_report(second_sites)

        # Align two motifs (using dynamic prgramming)
        alignment = motif_alignment(first_ensemble_report.aligned_sites,
                                    second_ensemble_report.aligned_sites)
        aligned_weblogos = (bioutils.weblogo_uri(alignment[0]),
                            bioutils.weblogo_uri(alignment[1]))

        # Clear session data from previous requests
        request.session['permuted_motifs'] = None
        request.session.modified = True

        return render(
            request, 'compare_motifs_results.html',
            {'form_title': "Step 3 of 3: Motif comparison results",
             'first_motif_reports': first_motif_reports,
             'first_ensemble_report': first_ensemble_report,
             'second_motif_reports': second_motif_reports,
             'second_ensemble_report': second_ensemble_report,
             'aligned_weblogos': aligned_weblogos})


def motif_similarity_measure(request):
    """AJAX handler for motif similarity measurements.

    Request has the similarity function and two motifs (space separated
    strings). Returns HTML which will be embedded into the main comparison
    page.

    The motif similarity function can be
    - Euclidean distance
    - Pearson correlation coeficient
    - Average log-likelihood ratio
    - Kullback-Leibler divergence
    - Levenshtein (edit) distance
    """

    unaligned_sites_first = request.POST['first_unaligned'].strip().split()
    unaligned_sites_second = request.POST['second_unaligned'].strip().split()
    sites_first = request.POST['first_aligned'].strip().split()
    sites_second = request.POST['second_aligned'].strip().split()

    # Site-based similarity measure (i.e. Levenshtein distance)
    if request.POST['similarity_function'] == 'site_based':
        boxplot, histogram = levenshtein_measure(sites_first, sites_second)
        return render(request, 'motif_similarity_levenshtein.html',
                      {'boxplot': boxplot,
                       'histogram': histogram,
                       'first_unaligned': unaligned_sites_first,
                       'second_unaligned': unaligned_sites_second,
                       'first_aligned': sites_first,
                       'second_aligned': sites_second})

    else:
        # Other similarity metrics
        functions = {
            'euclidean': euclidean_distance,
            'pearson': pearson_correlation_coefficient,
            'kullback_leibler': kullback_leibler_divergence,
            'average_log_likelihood': average_log_likelihood_ratio}
        similarity_function = functions[request.POST['similarity_function']]
        # Generate the histogram of scores and calculate the p-value
        fig, p_value = permutation_test(request, sites_first, sites_second,
                                        similarity_function)
        return render(
            request,
            'motif_similarity_%s.html' % request.POST['similarity_function'],
            {'plot': fig,
             'score': '%.3lf' % similarity_function(sites_first, sites_second),
             'p_value': '%.3lf' % p_value})


def levenshtein_distances(motif_a, motif_b):
    """Levenshtein distances between all pairss of sites from two motifs.

    Given two sets of binding sites (motif_a and motif_b), measures Levenshtein
    distance between all pairs of sequences from two motifs, to measure the
    distance between two motifs.
    """
    def minimum_edit_distance(s1, s2):
        """Edit distance between two sequences."""
        if len(s1) > len(s2):
            s1, s2 = s2, s1
        distances = range(len(s1) + 1)
        for index2, char2 in enumerate(s2):
            new_distances = [index2+1]
            for index1, char1 in enumerate(s1):
                if char1 == char2:
                    new_distances.append(distances[index1])
                else:
                    new_distances.append(1 + min((distances[index1],
                                                  distances[index1+1],
                                                  new_distances[-1])))
            distances = new_distances
        return distances[-1]

    dists = [float(minimum_edit_distance(seqa, seqb)) / (len(seqa) * len(seqb))
             for seqa in motif_a for seqb in motif_b]
    return dists


def levenshtein_measure(motif_a, motif_b):
    """Levenshtein distance between two motifs."""

    a_vs_a = levenshtein_distances(motif_a, motif_a)
    a_vs_b = levenshtein_distances(motif_a, motif_b)
    b_vs_b = levenshtein_distances(motif_b, motif_b)

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
    histogram = fig2img(plt.gcf())
    return boxplot, histogram


def motif_count_columns(sites):
    """Returns the columns of the count matrix of the motif."""
    motif = bioutils.build_motif(sites)
    motif.pseudocounts = 1.0
    return [{let: motif.counts[let][i] for let in 'ACGT'}
            for i in xrange(len(motif))]


def motif_pwm_columns(sites):
    """Returns the columns of teh PWM of the motif."""
    motif = bioutils.build_motif(sites)
    motif.pseudocounts = 1.0
    return [{let: motif.pwm[let][i] for let in 'ACGT'}
            for i in xrange(len(motif))]


def pearson_correlation_coefficient(sites_a, sites_b):
    """Returns the Pearson correlation coefficient of two motifs."""
    def pcc(cola, colb):
        cola_avg = sum(cola[let] for let in 'ACTG') / 4.0
        colb_avg = sum(colb[let] for let in 'ACTG') / 4.0
        return (sum(((cola[let]-cola_avg) * (colb[let]-colb_avg))
                    for let in 'ACTG') /
                math.sqrt(sum((cola[let]-cola_avg)**2 for let in 'ACTG') *
                          sum((colb[let]-colb_avg)**2 for let in 'ACTG')))

    colsa = motif_pwm_columns(sites_a)
    colsb = motif_pwm_columns(sites_b)
    return sum(pcc(*cols) for cols in zip(colsa, colsb))


def average_log_likelihood_ratio(sites_a, sites_b):
    """Returns the Average Log-likelihood ratio distance of two motfis."""

    def allr(cola, colb, cnta, cntb):
        return (sum((cnta[let]*safe_log2(colb[let]/0.25) +
                     cntb[let]*safe_log2(cola[let]/0.25))
                    for let in 'ACTG') /
                sum(cnta[let] + cntb[let] for let in 'ACTG'))

    return sum(allr(*args) for args in zip(
        motif_pwm_columns(sites_a), motif_pwm_columns(sites_b),
        motif_count_columns(sites_a), motif_count_columns(sites_b)))


def euclidean_distance(sites_a, sites_b):
    """Returns Euclidean distance between two sets of sites"""
    def ed(cola, colb):
        return math.sqrt(sum((cola[let] - colb[let])**2 for let in "ACGT"))
    return sum(ed(*cols) for cols in zip(motif_pwm_columns(sites_a),
                                         motif_pwm_columns(sites_b)))


def kullback_leibler_divergence(sites_a, sites_b):
    """Returns Kullback-Leibler divergence between two sets of sites."""
    def kl(cola, colb):
        return (sum(cola[let] * safe_log2(cola[let] / colb[let])
                    for let in 'ACTG') +
                sum(colb[let] * safe_log2(colb[let] / cola[let])
                    for let in 'ACTG')) / 2.0

    return sum(kl(*cols) for cols in zip(motif_pwm_columns(sites_a),
                                         motif_pwm_columns(sites_b)))


def motif_alignment(sites_a, sites_b):
    """Aligns two motifs that gives the maximum similarity.

    Returns two sets of sites that are aligned and of same length.
    """
    return motif_SW(sites_a, sites_b)


def motif_SW(sites_a, sites_b):
    """Smith-Waterman alignment of two motifs."""

    len_motif_a = len(sites_a[0])
    len_motif_b = len(sites_b[0])
    M = [[None for i in xrange(len_motif_b+1)] for j in xrange(len_motif_a+1)]
    for j in xrange(len_motif_b+1):
        M[0][j] = 0
    for i in xrange(len_motif_a+1):
        M[i][0] = 0
    # Fill the matrix
    for i in xrange(1, len_motif_a+1):
        for j in xrange(1, len_motif_b+1):
            M[i][j] = max(0, (M[i-1][j-1] +
                              pearson_correlation_coefficient(
                                  [site[i-1] for site in sites_a],
                                  [site[j-1] for site in sites_b])))

    # find max/min
    opt_score = M[0][0]
    opt_i = 0
    opt_j = 0
    for i in xrange(len(M)):
        for j in xrange(len(M[0])):
            if M[i][j] > opt_score:
                opt_score = M[i][j]
                opt_i, opt_j = i, j

    alignment_a_start = opt_i-1
    alignment_a_end = opt_i
    alignment_b_start = opt_j-1
    alignment_b_end = opt_j

    while opt_i >= 0 and opt_j >= 0 and M[opt_i][opt_j] > 0:
        alignment_a_start = opt_i-1
        alignment_b_start = opt_j-1
        opt_i -= 1
        opt_j -= 1

    return ([site[alignment_a_start:alignment_a_end] for site in sites_a],
            [site[alignment_b_start:alignment_b_end] for site in sites_b])


def permutation_test(request, ma, mb, fnc):
    """Performs the permutation test for the given motifs."""
    permuted_dists = permuted_distances(request, ma, mb, fnc)
    true_dist = fnc(ma, mb)
    plt.hist(permuted_dists, bins=30, normed=False, color='#0088CC',
             label='permuted pairs')
    plt.axvline(true_dist, linestyle='dashed', linewidth=2, color='#FF3300',
                label='motif pair')
    plt.xlabel({
        euclidean_distance: 'Eucledian distance',
        pearson_correlation_coefficient: 'Pearson Correlation Coefficient',
        kullback_leibler_divergence: 'Kullback-Leibler Divergence',
        average_log_likelihood_ratio: 'average log likelihood ratio'}[fnc])
    plt.ylabel('frequency')
    plt.legend()
    # calc p-value
    if fnc in [euclidean_distance, kullback_leibler_divergence]:
        p_value = (
            (sum(1 if pd < true_dist else 0 for pd in permuted_dists) + 1.0) /
            (len(permuted_dists) + 1))
    else:
        p_value = (
            (sum(1 if pd > true_dist else 0 for pd in permuted_dists) + 1.0) /
            (len(permuted_dists)+1))
    return fig2img(plt.gcf()), p_value


def sample(n, xs, replace=True):
    """Samples n objects from the list xs."""
    if replace:
        return [random.choice(xs) for i in range(n)]
    else:
        ys = list(xs[:])
        samp = []
        for i in range(n):
            y = random.choice(ys)
            samp.append(y)
            ys.remove(y)
        return samp


def permute_sites(sites):
    """Permutes columns of the binding motif."""
    sites = sites[:]
    l = len(sites[0])
    p = sample(l, range(l), replace=False)
    for i, site in enumerate(sites):
        sites[i] = "".join(site[p[j]] for j in p)
    return sites


def permuted_distances(request, motif_a, motif_b, dist_fun, n=100):
    """Permutes columns of two motifs and measures the distances.

    Performs this for n times and returns the list of scores."""

    # Don't permute if already done.
    if not request.session.get('permuted_motifs', None):
        request.session['permuted_motifs'] = (
            [permute_sites(motif_a) for i in xrange(n)],
            [permute_sites(motif_b) for i in xrange(n)])
        request.session.modified = True
        print 'calcluating'
    dists = [dist_fun(a, b)
             for a, b in zip(*request.session['permuted_motifs'])]
    return dists


def fig2img(fig):
    """Returns data URI, given the matplotlib plot."""
    imgdata = StringIO.StringIO()
    fig.savefig(imgdata, format='png')
    plt.clf()
    imgdata.seek(0)             # rewind the data
    return "data:image/png;base64,%s" % b64encode(imgdata.buf)


def safe_log2(x):
    return math.log(x, 2) if x != 0 else 0.0
