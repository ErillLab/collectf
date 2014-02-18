from django.shortcuts import render
from django.shortcuts import render_to_response
from django.template import RequestContext


import models
import random
from browse import motif_report

def about(request):
    return render_to_response("about.html", {},
                              context_instance=RequestContext(request))

def browse(request):
    return render_to_response("browse.html", {},
                              context_instance=RequestContext(request))

def search(request):
    return render_to_response("search.html", {},
                              context_instance=RequestContext(request))

def compare(request):
    return render_to_response("compare.html", {},
                              context_instance=RequestContext(request))

def contribute(request):
    return render_to_response("contribute.html", {},
                              context_instance=RequestContext(request))

def feedback(request):
    return render_to_response("feedback.html", {},
                              context_instance=RequestContext(request))

def stats(request):
    return render_to_response("stats.html", {},
                              context_instance=RequestContext(request))

def cite(request):
    return render_to_response("cite.html", {},
                              context_instance=RequestContext(request))

def links(request):
    return render_to_response("links.html", {},
                              context_instance=RequestContext(request))

def acknowledgements(request):
    return render_to_response("acknowledgements.html", {},
                              context_instance=RequestContext(request))

def greet(request):
    """Handler for the main page!
    It grabs a random motif from the database to display on the main page."""
    random_record = get_random_motif()
    template_dict = dict(TF_name = random_record.TF_name,
                         TF_accession = random_record.TF_accession,
                         organism = random_record.species_name,
                         genome_accession = random_record.genome_accession,
                         aligned_sites = random_record.align_sites(),
                         cur_site_insts = random_record.get_all_cur_site_insts(),
                    )

    return render_to_response("greet.html", {'random_rec': template_dict},
                              context_instance=RequestContext(request))

def get_random_motif(motif_len_th=30, motif_sz_th=10):
    """Get random motif from the database, to display on the main page. Randomly
    selected motif must be no longer than <motif_len_th> and motif size must be
    at least <motif_sz_th>."""
    # Get all possible TF-instance & species combinations
    TF_genome_list = models.Curation_SiteInstance.objects.values_list(
            'curation__TF_instances',
            'site_instance__genome').distinct()
    # select one of the combinations that satisfies the criteria (motif_len and
    # motif_sz)
    while True:
        TF_instances,genome = random.choice(TF_genome_list)
        # get all curation-site-instance objects for the TF-instance and genome
        cur_site_insts = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances = TF_instances,
            site_instance__genome = genome,
            site_type="motif_associated")
        # generate a motif report out of curation-site-instance objects
        if cur_site_insts.count() > motif_sz_th:
            # First criterion is satisfied.
            report = motif_report.MotifReport(cur_site_insts)
            # Align binding sites
            aligned_sites = report.align_sites()
            if len(aligned_sites[0]) < motif_len_th:
                # Second criterion is satisfied too, incorporate
                # non-motif-associated data into the motif-report.
                non_motif_cur_site_insts = models.Curation_SiteInstance.objects.filter(
                    curation__TF_instances = TF_instances,
                    site_instance__genome = genome,
                    site_type = "non_motif_associated")
                report.set_non_motif_curation_site_instances(non_motif_cur_site_insts)

                # An appropriate motif is randomly selected at this point. Return it.
                return report
            
