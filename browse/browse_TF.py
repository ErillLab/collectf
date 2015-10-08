"""This file contains the view functions for browsing by TF"""

from django.shortcuts import render_to_response
from django.template import RequestContext
from django.shortcuts import get_object_or_404, get_list_or_404
import browse.models as models
import browse.motif_report as motif_report

def browse_TF(request):
    """View function for browse by TF. Returns the TF tree"""
    TF_families = models.TFFamily.objects.all().order_by('name')
    return render_to_response('browse_TF.html', {'TF_families': TF_families},
                              context_instance=RequestContext(request))

def get_TFs(tf_type, tf_id):
    """Given a TF or family, return that TF or all TFs in that family"""
    tf_objs = models.TF.objects
    tfs = None
    if tf_type == 'all':
        tfs = tf_objs.all()
    elif tf_type == 'family':
        tf_family = get_object_or_404(models.TFFamily, TF_family_id=tf_id)
        tfs = get_list_or_404(models.TF, family=tf_family)
    elif tf_type == 'tf':
        tfs = get_list_or_404(models.TF, TF_id=tf_id)
    return tfs

def get_tf_instance(accession):
    return get_object_or_404(models.TFInstance, protein_accession=accession)

def get_results_TF(request, type_, id_):
    """GIven the type (TF or TF family) and the id of the object, return query
    results and list TF/species that have binding site data for the selected TF
    or TF family."""

    TFs = get_TFs(type_, id_)

    if type_ == 'tf':
        name = TFs[0].name
        desc = TFs[0].description
    elif type_ == 'family':
        # Get all TFs under that family
        TF_family = get_object_or_404(models.TFFamily, TF_family_id=id_)
        name = TF_family.name
        desc = TF_family.description

    # get all curation-site-instance objects for browsed TFs
    assert TFs
    cur_site_insts = models.Curation_SiteInstance.objects.\
                     filter(curation__TF_instances__TF__in=TFs)

    # generate all reports
    reports = motif_report.make_reports(cur_site_insts)

    return render_to_response("browse_results.html",
                              dict(title=name,
                                   description=desc,
                                   reports=[report.generate_browse_result_dict()
                                            for report in reports],
                                   tax_param_type='all',
                                   tax_param=-1,
                                   tf_param_type=type_,
                                   tf_param=id_,
                                   tech_param_type='all',
                                   tech_param=-1),
                              context_instance=RequestContext(request))

