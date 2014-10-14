from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import user_passes_test
from django.core.urlresolvers import reverse
from django.contrib import messages
from base import models
from base import bioutils
from django.utils.safestring import mark_safe

@user_passes_test(lambda u: u.is_superuser)
def home(request):
    """Entry point for validate curation"""
    template_file = "validate_curation_main.html"
    curations = models.Curation.objects.order_by('-created').all()
    return render_to_response(template_file, {'curations': curations},
                              context_instance=RequestContext(request))

@user_passes_test(lambda u: u.is_superuser)
def validate_curation(request, curation_id):
    """Validate curation handler."""
    if request.method == 'POST':
        curation = models.Curation.objects.get(pk=curation_id)
        curator = models.Curator.objects.get(user=request.user)
        curation.validated_by = curator
        curation.save()
        messages.add_message(request, messages.SUCCESS,
                             "Curation was validated successfully.")
        return HttpResponseRedirect(reverse(home))
    else:
        curation = models.Curation.objects.get(pk=curation_id)
        # Get motifs for this TF/species
        genomes = set(cs.site_instance.genome
                      for cs in curation.curation_siteinstance_set.all())
        csis = models.Curation_SiteInstance.objects.filter(
            site_type='motif_associated',
            site_instance__genome__in=genomes,
            curation__TF_instances=curation.TF_instances.all())
        motif_ids = csis.values_list('motif_id', flat=True).distinct()
        all_motifs = []
        for motif_id in motif_ids:
            sites = [csi.site_instance
                     for csi in csis.filter(motif_id=motif_id)]
            weblogo = mark_safe("<img src='%s'" %
                                bioutils.weblogo_uri(bioutils.run_lasagna(sites)))
            all_motifs.append((motif_id, weblogo))
        print len(all_motifs)
        # All motifs for this TF/species
        # initial data preparation functions
    return render(request, "validate_curation.html",
                  {'curation': curation, 'all_motifs': all_motifs},
                  context_instance=RequestContext(request))

