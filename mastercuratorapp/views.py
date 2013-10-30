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
from django.contrib.auth.decorators import user_passes_test


@user_passes_test(lambda u: u.is_superuser)
def home(request):
    template_file = "main.html"
    # get all not validated curations
    curations = models.Curation.objects.\
                filter(master_curator_verified=False).\
                order_by('created').all()
    template = {'curations': curations}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))

@user_passes_test(lambda u: u.is_superuser)
def validate_curation(request, curation_id):
    if request.method == 'POST':
        curation = models.Curation.objects.get(pk=request.POST["curation_id"])
        form = forms.EditCurationForm(request.POST, curation=curation)
        
        if form.is_valid():
            form_done(form, curation)
            messages.add_message(request, messages.SUCCESS, "Curation was modified/validated successfully.")
            return HttpResponseRedirect(reverse(browseapp.view_curation.view_curation, kwargs={'cid': curation.pk}))
    else:
        curation = models.Curation.objects.get(pk=curation_id)
        # initial data preparation functions
        form = get_form(curation)
        template = {'form': form,
                    'curation': curation}
        
    return render(request, "validate_curation.html",
                  {'form': form, 'curation': curation},
                  context_instance=RequestContext(request))

def get_form(curation):
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

def form_done(form, curation):
    """Called when the form is submitted."""
    cd = form.cleaned_data
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

    def site_instances_done():
        for k,v in cd.items():
            if k.startswith('site_instance'):
                site_id = int(k.replace("site_instance_", ""))
                curation_site_instance = models.Curation_SiteInstance.objects.get(pk=site_id)
                if v == False:
                    curation_site_instance.delete()

    TF_genome_done()
    techniques_done()
    curation_review_done()
    site_instances_done()
    curation.save()




