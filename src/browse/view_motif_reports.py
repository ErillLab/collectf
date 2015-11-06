"""This module contains function that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
Curation_SiteInstance objects, they are grouped by TF and species and rendered
properly.
"""

from django.shortcuts import get_object_or_404
from django.shortcuts import render

from core import models


def group_curation_site_instances(curation_site_instances):
    """Groups Curation_SiteInstance objects by TF-instance and genome."""
    TF_genome_pairs = curation_site_instances.values_list(
        'site_instance__genome', 'curation__TF_instances').distinct()
    return [curation_site_instances.filter(site_instance__genome=genome,
                                           curation__TF_instances=TF_instances)
            for genome, TF_instances in TF_genome_pairs]


def view_reports_by_TF_family(request, object_id):
    """Returns the motif reports for the given TF family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    curation_site_instance_grps = group_curation_site_instances(
        models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF__family=TF_family))
    return render(request, 'view_motif_reports.html',
                  {'csi': curation_site_instance_grps})
