from django.shortcuts import render_to_response
from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import user_passes_test
from django.core.urlresolvers import reverse
from django.contrib import messages
import forms
import django.forms as dforms
from base import models
from base.templatetags import diagram
from browse import view_curation

@user_passes_test(lambda u: u.is_superuser)
def home(request):
    """Entry point for validate curation"""
    template_file = "validate_curation_main.html"
    return render_to_response(template_file, {}, context_instance=RequestContext(request))

def view_edit_validate(request):
    """View main page for edit/validate curations"""
    template_file = "validate_edit_or_validate.html"
    # get all not validated curations
    curations = models.Curation.objects.\
                filter(validated_by=None).\
                order_by('created').all()
    template = {'curations': curations}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))

def view_validated_curations(request):
    """View page to view validated curations."""
    template_file = "validate_view_curations.html"
    # get all validated curation objects
    curations = models.Curation.objects.\
                exclude(validated_by=None).\
                order_by('created').all()
    template = {'curations': curations}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))

@user_passes_test(lambda u: u.is_superuser)
def edit_curation(request, curation_id):
    return edit_curation.edit_curation(request, curation_id)


@user_passes_test(lambda u: u.is_superuser)
def validate_curation(request, curation_id):
    """Validate curation handler."""
    if request.method == 'POST':
        curation = models.Curation.objects.get(pk=request.POST["curation_id"])
        form = forms.EditCurationForm(request.POST, curation=curation)
        
        if form.is_valid():
            validate_form_done(request, form, curation)
            messages.add_message(request, messages.SUCCESS, "Curation was modified/validated successfully.")
            return HttpResponseRedirect(reverse(browseapp.view_curation.view_curation, kwargs={'cid': curation.pk}))
    else:
        curation = models.Curation.objects.get(pk=curation_id)
        # initial data preparation functions
        form = validate_get_form(curation)
        template = {'form': form,
                    'curation': curation}
        
    return render(request, "validate_curation.html",
                  {'form': form, 'curation': curation},
                  context_instance=RequestContext(request))

@user_passes_test(lambda u: u.is_superuser)
def edit_validated_curation(request, curation_id):
    if request.method == 'POST':
        curation = models.Curation.objects.get(pk=request.POST["curation_id"])
        form = forms.EditCurationForm(request.POST, curation=curation)
        if form.is_valid():
            validate_form_done(request, form, curation)
            messages.add_message(request, messages.SUCCESS, "Curation was modified/validated successfully.")
            return HttpResponseRedirect(reverse(browseapp.view_curation.view_curation, kwargs={'cid': curation.pk}))
    else:
        curation = models.Curation.objects.get(pk=curation_id)
        form = validate_get_form(curation)
        # check if any of the sites has been submitted to NCBI.
        ncbi_submitted = False
        csis = curation.curation_siteinstance_set.all()
        if models.NCBISubmission.objects.filter(curation_site_instance__in=csis):
            ncbi_submitted = True

        template = {'form': form,
                    'curation': curation,
                    'ncbi_submitted': ncbi_submitted,}
        
    return render(request, "validate_curation.html",
                  template,
                  context_instance=RequestContext(request))


def validate_get_form(curation):
    def get_genome_accession():
        if curation.site_instances.all():
            return curation.site_instances.all()[0].genome.genome_accession
        return ""
    def get_used_techniques():
        ts = curation.experimental_techniques.all()
        return [str(t.technique_id) for t in ts]
    def get_external_db():
        try:
            external_db = models.Curation_ExternalDatabase.objects.get(curation=curation)
        except models.Curation_ExternalDatabase.DoesNotExist:
            external_db = None
        return external_db

    def populate_site_instances(form):
        for csi in curation.curation_siteinstance_set.all():
            site_instance = csi.site_instance
            seq = csi.site_instance.seq
            strand = '+' if site_instance.strand == 1 else '-'
            loc = '[%d,%d]' % (site_instance.start+1, site_instance.end+1)
            label = pretty_print.site2label(csi.pk, seq+' '+strand+loc)
            help_text = gene_diagram.regulation_diagram(csi.regulation_set.all(), csi.site_instance)
            form.fields["site_instance_%d"%csi.pk] = dforms.BooleanField(label=label,
                                                                         help_text=help_text,
                                                                         required=False)
            form.fields["site_instance_%d"%csi.pk].initial = True
            
    external_db = get_external_db()
    data = dict(
        # genome/TF initialization
        TF = curation.TF,
        TF_type = curation.TF_type,
        genome_accession = get_genome_accession(),
        TF_accession = curation.TF_instance.protein_accession,
        TF_species = curation.TF_species,
        site_species = curation.site_species,
        # techniques initialization
        techniques = get_used_techniques(),
        experimental_process = curation.experimental_process,
        forms_complex = curation.forms_complex,
        complex_notes = curation.complex_notes,
        external_db_type = (external_db.external_database.ext_database_id
                            if external_db else None),
        external_db_accession = (external_db.accession_number
                                 if external_db else ""),
        # curation review initialization
        revision_reasons = curation.requires_revision,
        confidence = curation.confidence,
        paper_complete = curation.publication.curation_complete,
        NCBI_submission_ready = curation.NCBI_submission_ready,
        notes = curation.notes,
    )
    #form = forms.EditCurationForm(data)
    # add sites
    #populate_site_instances(form)
    kwargs = {'curation': curation}
    form = forms.EditCurationForm(data, **kwargs)
    return form
