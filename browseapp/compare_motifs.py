from browse_base import *
from django import forms
from django.core.exceptions import ValidationError
from django.contrib.formtools.wizard.views import CookieWizardView
import search
import view_results
import StringIO
import matplotlib.pyplot as plt
from matplotlib import rc
from base64 import b64encode

FORM_TITLES = ["Step 1/3: Search for the first motif.",
               "Step 2/3: Search for the second motif.",
               "Step 3/3: Motif comparison results"]

FORM_DESCRIPTIONS = [
"""Search in CollecTF is fully customizable. Just select a taxonomic unit (e.g. the
Vibrio genus), a transcription factor family or instance (e.g. LexA) and the set of
experimental techniques that reported sites should be backed by and proceed.""",

"""Search for the second motif to be compared. Just select a taxonomic unit (e.g. the
Vibrio genus), a transcription factor family or instance (e.g. LexA) and the set of
experimental techniques that reported sites should be backed by and proceed.""",

"",
]

class MotifComparisonForm(forms.Form):
    pass

FORMS = [("motif_a", MotifComparisonForm),
         ("motif_b", MotifComparisonForm),
         #("search_results", MotifComparisonForm),
         ("comparison_results", MotifComparisonForm)]

TEMPLATES = {"motif_a": "motif_compare_search.html",
             "motif_b": "motif_compare_search.html",
             "comparison_results": "motif_compare_comparison_results.html"}

def motif_comparison_get(request, step):
    template = search.search_get_template()
    template["form_title"] = FORM_TITLES[step-1]
    template["form_description"] = FORM_DESCRIPTIONS[step-1]
    return render_to_response("motif_compare_search.html",
                              template,
                              context_instance=RequestContext(request))

def motif_comparison_post(request, step):
    try:
        motif_csis, non_motif_csis = search.search_post_helper(request)
    except:
        message = "Please select at least one TF, species and experimental technique to search database."
        messages.add_message(request, messages.ERROR, message)
        return HttpResponseRedirect(reverse(eval('motif_comparison_step%d' % step)))
    if not motif_csis:
        message = "No results found for selected TF and species"
        messages.add_message(request, messages.ERROR, message)
        return HttpResponseRedirect(reverse(eval('motif_comparison_step%d' % step)))

    search_results = search.group_search_results(motif_csis, non_motif_csis)
    # store search results
    request.session['cstep%d' % step] = search_results
    request.session.modified = True

    if step == 1:
        return HttpResponseRedirect(reverse(motif_comparison_step2))
    else: # step2
        # render results
        # get motif_associated and non-motif_associated curation site instances
        motif_a_csi_list = request.session["cstep1"]["view_all_csis"]
        motif_a_ncsi_list = request.session["cstep1"]["view_all_ncsis"]
        motif_b_csi_list = request.session["cstep2"]["view_all_csis"]
        motif_b_ncsi_list = request.session["cstep2"]["view_all_ncsis"]
        # get the ensemble views
        motif_a_data = view_results.prepare_results(motif_a_csi_list, motif_a_ncsi_list)
        motif_b_data = view_results.prepare_results(motif_b_csi_list, motif_b_ncsi_list)
        sites_a = motif_a_data['ensemble_report']['aligned_sites']
        sites_b = motif_b_data['ensemble_report']['aligned_sites']
        aligned_sites_a, aligned_sites_b = motif_alignment(sites_a, sites_b)
        # clean session data from possible previous requests
        request.session['permute_motif_a'] = None
        request.session['permute_motif_b'] = None
        request.session.modified = True
        

        return render_to_response("motif_compare_comparison_results.html",
                                  {"form_title": FORM_TITLES[2],
                                   "form_description": FORM_DESCRIPTIONS[2],
                                   "motif_a": motif_a_data,
                                   "motif_b": motif_b_data,
                                   "motif_a_reports": request.session['cstep1']['reports'],
                                   "motif_b_reports": request.session['cstep2']['reports'],
                                   "aligned_sites_a": aligned_sites_a,
                                   "aligned_sites_b": aligned_sites_b},
                                  context_instance=RequestContext(request))
        
def motif_comparison_step1(request):
    if request.POST:
        return motif_comparison_post(request, 1)
    return motif_comparison_get(request, 1)

def motif_comparison_step2(request):
    if request.POST:
        return motif_comparison_post(request, 2)
    return motif_comparison_get(request, 2)

def motif_sim_measure(request):
    """
    AJAX handler for motif similarity measurements.
    Request has the similarity function and two motifs (comma seperated strings).
    Reponse returns HTML which is embedded into the main comparison page.
    """
    unaligned_sites_a = request.POST['unaligned_sites_a'].strip().split(',')
    unaligned_sites_b = request.POST['unaligned_sites_b'].strip().split(',')
    sites_a = request.POST['sites_a'].strip().split(',')
    sites_b = request.POST['sites_b'].strip().split(',')

    if request.POST['fun'] == 'site_based':
        boxplot, hist = levenshtein_measure(sites_a, sites_b)
        return render_to_response("motif_sim_levenshtein.html",
                                  {'boxplot': boxplot,
                                   'hist': hist,
                                   'unaligned_sites_a': unaligned_sites_a,
                                   'unaligned_sites_b': unaligned_sites_b,
                                   'aligned_sites_a': sites_a,
                                   'aligned_sites_b': sites_b},
                                  context_instance=RequestContext(request))
        
    else: # other similarity metrics
        try:
            fun_str = request.POST['fun']
            fun2call = {'ED': euclidean_distance,
                        'PCC': pearson_correlation_coefficient,
                        'KL': kullback_leibler_divergence,
                        'ALLR': average_log_likelihood_ratio,
                        }[fun_str]
            
            sites_a, sites_b = motif_alignment(sites_a, sites_b) # align motif_a and motif_b
            fig, p_val = motif_sim_test(request, sites_a, sites_b, fun2call)
        except Exception as e:
            print e
        return render_to_response("motif_sim_%s.html" % fun_str,
                                  {'plot': fig,
                                   'sc': '%.3lf' % fun2call(sites_a, sites_b),
                                   'p_value': '%.3lf' % p_val,
                                   'unaligned_sites_a': unaligned_sites_a,
                                   'unaligned_sites_b': unaligned_sites_b,
                                   'aligned_sites_a': sites_a,
                                   'aligned_sites_b': sites_b},
                                  context_instance=RequestContext(request))

def levenshtein_measure(motif_a, motif_b):
    a_vs_a = levenshtein_motifs(motif_a, motif_a)
    a_vs_b = levenshtein_motifs(motif_a, motif_b)
    b_vs_b = levenshtein_motifs(motif_b, motif_b)

    bp = plt.boxplot([a_vs_a, a_vs_b, b_vs_b])
    plt.xticks(range(1,4), [r'$M_a$ vs $M_a$', r'$M_a$ vs $M_b$', r'$M_b$ vs $M_b$'])
    plt.setp(bp['boxes'], color='#0088CC')
    plt.setp(bp['whiskers'], color='#0088CC')
    plt.setp(bp['fliers'], color='#0088CC', marker='+')
    plt.ylabel('Levenshtein distance')
    boxplot = fig2img(plt.gcf())
    plt.hist([a_vs_a, a_vs_b, b_vs_b], bins=15,
             label=[r'$M_a$ vs $M_a$', r'$M_a$ vs $M_b$', r'$M_b$ vs $M_b$'],
             color=["#FFCC33", "#006699", "#FF6633"])
    plt.ylabel('frequency of site-pairs')
    plt.xlabel('Levenshtein distance')
    plt.legend()
    hist = fig2img(plt.gcf())

    return boxplot, hist

def motif_sim_test(request, ma, mb, fnc):
    """Given two motifs and a similarity function, perform the permutation tests and
    return the histogram"""
    permuted_dists = permutation_test(request, ma, mb, fnc)
    true_dist = fnc(ma, mb)
    plt.hist(permuted_dists, bins=30, normed=False, color='#0088CC', label="permuted pairs")
    plt.axvline(true_dist, linestyle='dashed', linewidth=2, color='#FF3300', label="motif pair")
    plt.xlabel({euclidean_distance: "Eucledian distance",
                pearson_correlation_coefficient: "Pearson Correlation Coefficient",
                kullback_leibler_divergence: "Kullback-Leibler Divergence",
                average_log_likelihood_ratio: "average log likelihood ratio"}[fnc])
    plt.ylabel('frequency')
    plt.legend()
    # calc p-value
    if fnc in [euclidean_distance, kullback_leibler_divergence]:
        p_value = (sum(1 if pd<true_dist else 0 for pd in permuted_dists) + 1.0) / (len(permuted_dists)+1)
    else:
        p_value = (sum(1 if pd>true_dist else 0 for pd in permuted_dists) + 1.0) / (len(permuted_dists)+1)
    return fig2img(plt.gcf()), p_value

def motif_alignment(sites_a, sites_b):
    """
    Given two sets of sites and a similiarity/distance function, align two motifs
    that gives the maximum similarity (minimum distance). Returns two sets of sites
    that are aligned and of same length.
    """
    return motif_SW(sites_a, sites_b)
    
def motif_SW(sites_a, sites_b):
    # motif smith-waterman similarity
    len_motif_a = len(sites_a[0])
    len_motif_b = len(sites_b[0])
    M = [[None for i in xrange(len_motif_b+1)] for j in xrange(len_motif_a+1)]
    for j in xrange(len_motif_b+1):
        M[0][j] = 0
    for i in xrange(len_motif_a+1):
        M[i][0] = 0
    # fill the matrix
    for i in xrange(1, len_motif_a+1):
        for j in xrange(1, len_motif_b+1):
            # score from diag
            # ungapped
            M[i][j] = max(0,
                          (M[i-1][j-1] + pearson_correlation_coefficient([site[i-1] for site in sites_a],\
                                                                         [site[j-1] for site in sites_b])))

    # find max/min
    opt_score = M[0][0]
    opt_i = 0
    opt_j = 0
    for i in xrange(len(M)):
        for j in xrange(len(M[0])):
            if M[i][j] > opt_score:
                opt_score = M[i][j]
                opt_i, opt_j = i,j

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
    
    
def fig2img(fig):
    """Given matplotlib plot return data URI."""
    imgdata = StringIO.StringIO()
    fig.savefig(imgdata, format='png')
    plt.clf()
    imgdata.seek(0) # rewind the data
    return "data:image/png;base64,%s" % b64encode(imgdata.buf)

def sample(n,xs,replace=True):
    """Sample n objects from the list xs"""
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

def permute_motif(motif):
    """Permute columns of the binding motif"""
    motif = motif[:]
    l = len(motif[0])
    p = sample(l, range(l), replace=False)
    for i,site in enumerate(motif):
        motif[i] = "".join(site[p[j]] for j in p)
    return motif

def permutation_test(request, motif_a, motif_b, dist_fun, n=100):
    """Permute columns of two motifs and measure the similarity/distance with the
    specified function. Do this for n times and return the list of scores."""
    if 'permuted_motif_a' not in request.session:
        print 'creating new permutations for a'
        request.session['permuted_motif_a'] = [permute_motif(motif_a) for i in xrange(n)]
        request.session.modified = True
    if 'permuted_motif_b' not in request.session:
        print 'creating new permutations for b'
        request.session['permuted_motif_b'] = [permute_motif(motif_b) for i in xrange(n)]
        request.session.modified = True
    dists = [dist_fun(a,b) for a,b in zip(request.session['permuted_motif_a'],
                                          request.session['permuted_motif_b'])]
    return dists

def levenshtein(seq1, seq2):
    """Levenshtein distance between two sequences."""
    oneago = None
    thisrow = range(1, len(seq2) + 1) + [0]
    for x in xrange(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]

def levenshtein_motifs(motif_a, motif_b):
    """
    Given two sets of binding sites (motif_a and motif_b), measure Levenshtein
    distance between all pairs of sequences from two motifs, to measure the distance
    between two motifs
    """
    dists = [float(levenshtein(seqa, seqb)) / (len(seqa) * len(seqb))
            for seqa in motif_a for seqb in motif_b]
    return dists

def euclidean_distance(sites_a, sites_b):
    """Euclidean distance between two sets of sites"""
    ma = bioutils.create_motif(sites_a)
    mb = bioutils.create_motif(sites_b)
    def ed(cola, colb):
        return math.sqrt(sum((cola[l] - colb[l])**2 for l in "ACGT"))
    return sum(ed(cola, colb) for (cola,colb) in zip(ma.pwm(), mb.pwm()))

def kullback_leibler_divergence(sites_a, sites_b):
    """Kullback-Leibler divergence between two sets of sites"""
    def safe_log2(x):
        return math.log(x,2) if x != 0 else 0.0

    def kl(cola, colb):
        return (sum(cola[l] * safe_log2(cola[l] / colb[l]) for l in "ACTG") +
                sum(colb[l] * safe_log2(colb[l] / cola[l]) for l in "ACTG")) / 2.0

    ma = bioutils.create_motif(sites_a)
    mb = bioutils.create_motif(sites_b)
    return sum(kl(cola, colb) for (cola,colb) in zip(ma.pwm(), mb.pwm()))

def pearson_correlation_coefficient(sites_a, sites_b):
    """PEarson correlation coefficient"""
    def pcc(cola, colb):
        cola_avg = sum(cola[l] for l in "ACTG") / 4.0
        colb_avg = sum(colb[l] for l in "ACTG") / 4.0
        return (sum(((cola[l]-cola_avg) * (colb[l]-colb_avg)) for l in "ACTG") /
                math.sqrt(sum((cola[l]-cola_avg)**2 for l in "ACTG") *
                          sum((colb[l]-colb_avg)**2 for l in "ACTG")))

    ma = bioutils.create_motif(sites_a)
    mb = bioutils.create_motif(sites_b)
    return sum(pcc(cola, colb) for (cola,colb) in zip(ma.pwm(), mb.pwm()))

def average_log_likelihood_ratio(sites_a, sites_b):
    """Average Log-likelihood ratio distance"""
    def safe_log2(x):
        return math.log(x,2) if x != 0 else 0.0
    def allr(cola, colb, cnta, cntb):
        return (sum((cnta[l]*safe_log2(colb[l]/0.25) +
                     cntb[l]*safe_log2(cola[l]/0.25))
                    for l in "ACTG") /
                sum(cnta[l] + cntb[l] for l in 'ACTG'))
    ma = bioutils.create_motif(sites_a)
    mb = bioutils.create_motif(sites_b)
    # reformat biopython count matrices
    counts_a = [dict((l, ma.counts[l][i]) for l in "ACTG") for i in xrange(ma.length)]
    counts_b = [dict((l, mb.counts[l][i]) for l in "ACTG") for i in xrange(mb.length)]
    return sum(allr(cola, colb, cnta, cntb)
               for (cola,colb,cnta,cntb) in zip(ma.pwm(), mb.pwm(), counts_a, counts_b))
        
