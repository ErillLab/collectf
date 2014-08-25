"""This file contains function definitions that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
curation-site-instance objects, they are grouped by TF and species and rendered
properly."""

from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import Http404
import browse.models as models
import browse.motif_report as motif_report
import browse.browse_tax as browse_tax
import browse.browse_TF as browse_TF
import browse.browse_tech as browse_tech

def view_reports(request, tax_param_type, tax_param,
                 tf_param_type, tf_param,
                 tech_param_type, tech_param,
                 integrate_non_motif):

    """Given an organism and a TF, find all sites and make reports out of
    them.

    tax_param_type:  all/group/species
    tax_param:       id of the taxonomy object
    tf_param_type:   all/family/tf
    tf_param:        id of the TF/TF_family object
    tech_param_type: all/binding/expression/binding_category/expression_category
    integrate_non_motif: 1 to integrate non-motif-associated data, 0 otherwise
    """

    orgs = browse_tax.get_species(tax_param_type, tax_param)
    tfs = browse_TF.get_TFs(tf_param_type, tf_param)
    _, _, techs = browse_tech.get_techniques(tech_param_type, tech_param)

    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in=orgs,
        site_instance__curation__TF__in=tfs,
        experimental_techniques__in=techs,
    )

    # Exclude non-motif-associated sites here, it can be integrated if user
    # wants to do so
    if integrate_non_motif != 1:
        cur_site_insts = \
            cur_site_insts.exclude(site_type='non_motif_associated')

    reports = motif_report.make_distinct_reports(cur_site_insts)
    # Combine reports
    ensemble_report = motif_report.make_ensemble_report(cur_site_insts)

    if integrate_non_motif not in ['0', '1']:
        raise Http404
    switched_integrate = 1 if integrate_non_motif == '0' else 0

    return render_to_response('view_reports.html',
                              dict(reports=[report.generate_view_reports_dict()
                                            for report in reports],
                                   ensemble_report=
                                   ensemble_report.generate_view_reports_dict(),
                                   tax_param_type=tax_param_type,
                                   tax_param=tax_param,
                                   tf_param_type=tf_param_type,
                                   tf_param=tf_param,
                                   tech_param_type=tech_param_type,
                                   tech_param=tech_param,
                                   integrate_non_motif=switched_integrate),
                              context_instance=RequestContext(request))

def view_reports_by_id_list(request):
    """Given a collection of curation-site-instance ids with a post request,
    group them by TF/species and return the reports"""
    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        pk__in=request.POST['csi_list'].strip().split(',')
    )

    integrate_non_motif = bool('integrate_non_motif' in request.POST)

    if not integrate_non_motif:
        cur_site_insts =\
            cur_site_insts.exclude(site_type='non_motif_associated')

    reports = motif_report.make_distinct_reports(cur_site_insts)
    ensemble_report = motif_report.make_ensemble_report(cur_site_insts)

    return render_to_response('view_reports_by_id_list.html',
                              dict(reports=[report.generate_view_reports_dict()
                                            for report in reports],
                                   ensemble_report=
                                   ensemble_report.generate_view_reports_dict(),
                                   cur_site_insts=','.join(
                                       map(lambda csi: '%d'%csi.pk,
                                           cur_site_insts)),
                                   integrate_non_motif=integrate_non_motif),
                              context_instance=RequestContext(request))
