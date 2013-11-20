from django.shortcuts import render_to_response
from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages
import models
import forms
import django.forms as dforms
from baseapp.templatetags import pretty_print, gene_diagram
import browseapp.view_curation
import collectfapp.editcurationview
from django.contrib.auth.decorators import user_passes_test

@user_passes_test(lambda u: u.is_superuser)
def home(request):
    template_file = "main.html"
    return render_to_response(template_file, {}, context_instance=RequestContext(request))

def validate_edit_main(request):
    template_file = "validate_edit_main.html"
    # get all not validated curations
    curations = models.Curation.objects.\
                filter(validated_by=None).\
                order_by('created').all()
    template = {'curations': curations}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))

def view_validated_curations(request):
    template_file = "view_validated_curations.html"
    # get all validated curation objects
    curations = models.Curation.objects.\
                exclude(validated_by=None).\
                order_by('created').all()
    template = {'curations': curations}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))

@user_passes_test(lambda u: u.is_superuser)
def edit_curation(request, curation_id):
    return collectfapp.editcurationview.edit_curation(request, curation_id)

@user_passes_test(lambda u: u.is_superuser)
def validate_curation(request, curation_id):
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
        print csis
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
        TF_function = curation.TF_function,
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

def validate_form_done(request, form, curation):
    """Called when the form is submitted."""
    cd = form.cleaned_data
    assert request.method == "POST"
    assert "ncbi_submitted" in request.POST
    assert "why_obsolete" in request.POST
    cd['ncbi_submitted'] = True if request.POST['ncbi_submitted']=='True' else False
    cd['why_obsolete'] = request.POST['why_obsolete']
    print 'ncbi_submitted', cd['ncbi_submitted']
    print 'why_obsolete', cd['why_obsolete']
    def TF_genome_done():
        # TF/Genome update
        curation.TF = cd['TF']
        curation.TF_type = cd['TF_type']
        curation.TF_function = cd['TF_function']
        curation.TF_species = cd['TF_species']
        curation.site_species = cd['site_species']
        # TF_instance
        tf_ins = models.TFInstance.objects.get(protein_accession=cd['TF_accession'])
        curation.TF_instance = tf_ins

    def techniques_done():
        # Experimental techniques update
        curation.experimental_techniques.clear()
        for t in cd['techniques']:
            curation.experimental_techniques.add(t)
        curation.experimental_process = cd['experimental_process']
        curation.forms_complex = cd['forms_complex']
        curation.complex_notes = cd['complex_notes']
        # external db_type
        # remove existings
        models.Curation_ExternalDatabase.objects.filter(curation=curation).delete()
        # add new one
        if cd['external_db_type'] != "None" and cd['external_db_accession']:
            external_db_type = models.ExternalDatabase.objects.get(ext_database_id=cd['external_db_type'])
            curation_ext_ref = models.Curation_ExternalDatabase(curation=curation,
                                                                external_database=external_db_type,
                                                                accession_number=cd['external_db_accession'])
            curation_ext_ref.save()
            
    def curation_review_done():
        curation.requires_revision = cd['revision_reasons']
        curation.confidence = cd['confidence']
        curation.notes = cd['notes']
        curation.NCBI_submission_ready = cd['NCBI_submission_ready']
        publication = curation.publication
        if cd['paper_complete']:
            publication.curation_complete = True
        else:
            publication.curation_complete = False
        publication.save()

    def site_instances_done(ncbi_submitted, why_obsolete):
        # if a site is deselected
        # (a) if ncbi_submitted, don't delete it, mark it as obsolete
        # (b) otherwise, delete
        for k,v in cd.items():
            if k.startswith('site_instance'):
                site_id = int(k.replace("site_instance_", ""))
                curation_site_instance = models.Curation_SiteInstance.objects.get(pk=site_id)
                if v == False:
                    if ncbi_submitted:
                        curation_site_instance.is_obsolete = True
                        curation_site_instance.why_obsolete = why_obsolete
                        curation_site_instance.save()
                    else:
                        curation_site_instance.delete()

    TF_genome_done()
    techniques_done()
    curation_review_done()
    site_instances_done(cd['ncbi_submitted'], cd['why_obsolete'])
    curator = models.Curator.objects.get(user=request.user)
    print curator
    curation.validated_by = curator
    curation.save()
    
@user_passes_test(lambda u: u.is_superuser)
def withdraw_site(request):
    """Mark a site instance as obsolete"""
    assert request.POST
    from baseapp.templatetags import dbxref_utils
    site_id = request.POST['csi_id']
    why_obsolete = request.POST['why_obsolete']
    csi = models.Curation_SiteInstance.objects.get(pk=site_id)
    csi.is_obsolete = True
    csi.why_obsolete = why_obsolete
    csi.save()
    messages.add_message(request, messages.SUCCESS, "Site was marked as obsolete.")
    print dbxref_utils.id2dbxref(site_id)
    return browseapp.view_site.browse_by_site(request, dbxref_utils.id2dbxref_hex_only(site_id))
    
