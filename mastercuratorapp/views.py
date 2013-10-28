from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
import models
import forms
import django.forms as dforms
from baseapp.templatetags import pretty_print, gene_diagram

def home(request):
    template_file = "main.html"
    # get all not validated curations
    curations = models.Curation.objects.\
                filter(master_curator_verified=False).\
                order_by('created').all()
    template = {'curations': curations}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))

def validate_curation(request, curation_id):
    curation = models.Curation.objects.get(pk=curation_id)
    template_file = "validate_curation.html"
    # initial data preparation functions
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
            form.fields["site_instance_%d"] = dforms.BooleanField(label=label,
                                                                  help_text=help_text)
                                                      
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
    form = forms.EditCurationForm(data)
    # add sites
    populate_site_instances(form)
    template = {'form': form,
                'curation': curation}
    return render_to_response(template_file, template,
                              context_instance=RequestContext(request))
