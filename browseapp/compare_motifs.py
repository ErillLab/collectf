from browse_base import *
from django import forms
from django.contrib.formtools.wizard.views import CookieWizardView
import search
import view_results
import StringIO
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from base64 import b64encode
import random
import math

class MotifComparisonForm(forms.Form):
    pass

FORMS = [("motif_a", MotifComparisonForm),
         ("motif_b", MotifComparisonForm),
         #("search_results", MotifComparisonForm),
         ("comparison_results", MotifComparisonForm)]

TEMPLATES = {"motif_a": "motif_compare_search.html",
             "motif_b": "motif_compare_search.html",
             "comparison_results": "motif_compare_comparison_results.html"}

FORM_TITLES = {"motif_a": "Search for the first motif.",
               "motif_b": "Search for the second motif.",
               "comparison_results": "Motif comparison results"}

FORM_DESCRIPTIONS = {
    "motif_a": """Search in CollecTF is fully customizable. Just select a taxonomic
               unit (e.g. the Vibrio genus), a transcription factor family or
               instance (e.g. LexA) and the set of experimental techniques that
               reported sites should be backed by and proceed.""",

    "motif_b": """Search for the second motif to be compared. Just select a taxonomic
               unit (e.g. the Vibrio genus), a transcription factor family or
               instance (e.g. LexA) and the set of experimental techniques that
               reported sites should be backed by and proceed.""",

    "comparison_results": "",
    }

class MotifComparisonWizard(CookieWizardView):
    def get_template_names(self):
        return [TEMPLATES[self.steps.current]]

    def get_context_data(self, form, **kwargs):
        """
        Update context to include search form in the first two steps of motif
        comparison.
        """
        c = super(MotifComparisonWizard, self).get_context_data(form=form, **kwargs)
        # Update title and description
        c["form_title"] = FORM_TITLES[self.steps.current]
        c["form_description"] = FORM_DESCRIPTIONS[self.steps.current]

        # If it is motif search step, generate the search form
        if self.steps.current in ["motif_a", "motif_b"]:
            template = search.search_get_template()
            c.update(template)

        if self.steps.current == "search_results":
            c.update({"motif_a": self.request.session["motif_a"],
                      "motif_b": self.request.session["motif_b"]})

        if self.steps.current == "comparison_results":
            # get motif_associated and non-motif_associated curation site instances
            motif_a_csi_list = self.request.session["motif_a"]["view_all_csis"]
            motif_a_ncsi_list = self.request.session["motif_a"]["view_all_ncsis"]
            motif_b_csi_list = self.request.session["motif_b"]["view_all_csis"]
            motif_b_ncsi_list = self.request.session["motif_b"]["view_all_ncsis"]
            # get the ensemble views
            motif_a_data = view_results.prepare_results(motif_a_csi_list, motif_a_ncsi_list)
            motif_b_data = view_results.prepare_results(motif_b_csi_list, motif_b_ncsi_list)
            # put list of TF and species names in template
            get_TF_name = lambda reports: list(set(map(lambda rep: rep["TF_name"], reports)))
            get_sp_name = lambda reports: list(set(map(lambda rep: rep["species_name"], reports)))
            motif_a_data["TFs"] = get_TF_name(self.request.session["motif_a"]["reports"])
            motif_b_data["TFs"] = get_TF_name(self.request.session["motif_b"]["reports"])
            motif_a_data["species"] = get_sp_name(self.request.session["motif_a"]["reports"])
            motif_b_data["species"] = get_sp_name(self.request.session["motif_b"]["reports"])
            c.update({"motif_a": motif_a_data,
                      "motif_b": motif_b_data,
                      })
        return c

    def process_step(self, form):
        """
        Process data after each step.
        Used to perform searches for both motifs before comparison step
        """
        if self.steps.current in ["motif_a", "motif_b"]: # if one of the motif search steps
            # get motif and non-motif -associated sites
            motif_csis, non_motif_csis = search.search_post_helper(self.request)
            search_results = search.group_search_results(motif_csis, non_motif_csis)
            # store search results
            self.request.session[self.steps.current] = search_results
            self.request.session.modified = True
            
        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        return render_to_response("success.html")

def motif_sim_measure(request):
    """
    AJAX handler for motif similarity measurements.
    Request has the similarity function and two motifs (comma seperated strings).
    Reponse returns HTML which is embedded into the main comparison page.
    """
    def levenshtein_measure(motif_a, motif_b):
        a_vs_a = levenshtein_motifs(motif_a, motif_a)
        a_vs_b = levenshtein_motifs(motif_a, motif_b)
        b_vs_b = levenshtein_motifs(motif_b, motif_b)
        plt.boxplot([a_vs_a, a_vs_b, b_vs_b])
        plt.xticks(range(1,4), [r'$M_a$ vs $M_a$', r'$M_a$ vs $M_b$', r'$M_b$ vs $M_b$'])
        boxplot = fig2img(plt.gcf())
        plt.hist([a_vs_a, a_vs_b, b_vs_b], bins=15, normed=1,
                 label=[r'$M_a$ vs $M_a$', r'$M_a$ vs $M_b$', r'$M_b$ vs $M_b$'])
        plt.legend()
        hist = fig2img(plt.gcf())
        return boxplot, hist

    sites_a = request.POST['sites_a'].strip().split(',')
    sites_b = request.POST['sites_b'].strip().split(',')
    # remove gaps for now
    sites_a = [site.replace('-', 'A') for site in sites_a]
    sites_b = [site.replace('-', 'A') for site in sites_b]
    if request.POST['fun'] == 'levenshtein':
        boxplot, hist = levenshtein_measure(sites_a, sites_b)
        return render_to_response("motif_sim_levenshtein.html",
                                  {'boxplot': boxplot,
                                   'hist': hist,},
                                  context_instance=RequestContext(request))
        
    else: # other similarity metrics
        fun_str = request.POST['fun']
        fun2call = {'ED': euclidean_distance,
                    'PCC': pearson_correlation_coefficient,
                    'KL': kullback_leibler_divergence,
                    'ALLR': average_log_likelihood_ratio,
                    }[fun_str]
        fig = motif_sim_test(sites_a, sites_b, fun2call)
        return render_to_response("motif_sim_%s.html" % fun_str,
                                  {'plot': fig,
                                   'sc': fun2call(sites_a, sites_b)},
                                  context_instance=RequestContext(request))

def motif_sim_test(ma, mb, fnc):
    """Given two motifs and a similarity function, perform the permutation tests and
    return the histogram"""
    permuted_dists = permutation_test(ma, mb, fnc)
    true_dist = fnc(ma, mb)
    plt.hist(permuted_dists, bins=30, normed=1, color='c')
    plt.axvline(true_dist, linestyle='dashed', linewidth=2, color='b')
    return fig2img(plt.gcf())

def motif_alignment(sites_a, sites_b, fnc):
    """
    Given two motifs and a similiarity function, align two motifs that gives the
    maximum similarity. Returns two offset values, one for ma and one for mb. It
    means that motifs are aligned starting from columns ma[offset_a] and
    mb[offset_b].
    """
    pass
    

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

def permutation_test(motif_a, motif_b, dist_fun, n=250):
    """Permute columns of two motifs and measure the similarity/distance with the
    specified function. Do this for n times and return the list of scores."""
    dists = [dist_fun(permute_motif(motif_a), permute_motif(motif_b))
             for i in xrange(n)]
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
        
