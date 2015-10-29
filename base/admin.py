from django import forms
from django.contrib import admin
from django.contrib import messages
from django.db.models import get_app
from django.db.models import get_models
from django.http import HttpResponseRedirect
from django.shortcuts import render_to_response
from django.template import RequestContext

from base.models import *

class GenomeAdmin(admin.ModelAdmin):
    list_display = ('genome_accession', 'organism')

class CurationSiteInstanceInline(admin.StackedInline):
    model = Curation_SiteInstance
    filter_horizontal = ('experimental_techniques',)
    extra = 0
    
class CurationAdmin(admin.ModelAdmin):
    def PMID(self, obj):
        return obj.publication.pmid
    filter_horizontal = ('TF_instances',)
    list_display = ('curation_id', 'site_species', 'PMID', 'TF',
                    'curator', 'validated_by', 'NCBI_submission_ready')
    list_filter = ('curator', 'NCBI_submission_ready', 'requires_revision')
    ordering = ('-curation_id',)
    inlines = (CurationSiteInstanceInline,)

class CuratorAdmin(admin.ModelAdmin):
    def username(self, obj):
        return obj.user.username
    def email(self, obj):
        return obj.user.email
    def name(self, obj):
        return obj.user.first_name + ' ' + obj.user.last_name
    list_display = ('username', 'email', 'name')

class SiteInstanceAdmin(admin.ModelAdmin):
    list_display = ('site_id',)
    ordering = ('-site_id', )

class Curation_SiteInstanceAdmin(admin.ModelAdmin):
    list_filter = ('site_type',)
    filter_horizontal = ("experimental_techniques",)
    ordering = ('-id',)

class NotAnnotatedSiteInstanceAdmin(admin.ModelAdmin):
    filter_horizontal = ('experimental_techniques',)
    ordering = ('-id',)
    
class BatchAssignForm(forms.Form):
    _selected_action = forms.CharField(widget=forms.MultipleHiddenInput)
    curator = forms.ModelChoiceField(Curator.objects)
    
def assign_papers(modeladmin, request, queryset):
    """Assign papers to selected curator"""
    # check if any selected curation is complete
    for obj in queryset:
        if obj.curation_complete:
            messages.error(request, "Some of the selected papers are already curated. "
                           "Curated papers can not be reassigned.")
            return

    # reassignment step
    form = None
    if 'apply' in request.POST:
        form = BatchAssignForm(request.POST)
        if form.is_valid():
            curator = form.cleaned_data['curator']
            for obj in queryset:
                obj.assigned_to = curator
                obj.save()
            messages.success(request, "Papers successfully assigned to %s." % curator.user.username)
            return HttpResponseRedirect(request.get_full_path())

    if not form:
        form = BatchAssignForm(initial={'_selected_action': request.POST.getlist(admin.ACTION_CHECKBOX_NAME)})

    return render_to_response('admin/batch_paper_assign.html',
                              {'papers': queryset, 'form': form},
                              context_instance = RequestContext(request))
assign_papers.short_description = "Assign selected papers to a curator"

class PublicationAdmin(admin.ModelAdmin):
    list_display = ('publication_id', 'pmid', 'title', 'reported_species', 'reported_TF', 'assigned_to')
    list_filter = ('assigned_to', 'curation_complete')
    search_fields = ('pmid', 'reported_TF', 'reported_species')
    actions = [assign_papers]
    

class ExperimentalTechniqueAdmin(admin.ModelAdmin):
    list_display = ('name', )
    filter_horizontal = ('categories',)

class GeneAdmin(admin.ModelAdmin):
    list_display = ('gene_accession', 'name')
    list_filter = ('genome__genome_accession',)

class GenomeAdmin(admin.ModelAdmin):
    list_display = ('genome_accession', 'organism')

class TFFamilyAdmin(admin.ModelAdmin):
    list_display = ('name',)
    ordering = ('name',)
    
class TFAdmin(admin.ModelAdmin):
    list_display = ('name', 'family')
    list_filter = ('family',)
    ordering = ('name',)

class TFInstance(admin.ModelAdmin):
    list_display = ('protein_accession', 'name', 'TF')
    list_filter = ('')
    ordering = ('name')

class NCBISubmissionAdmin(admin.ModelAdmin):
    list_display = ('pk', 'genome_submitted_to', 'submission_time')

    
def register_rest():
    """Register all models that have not been registered explicitly."""
    app = get_app("base")
    for model in get_models(app):
        try:
            admin.site.register(model)
        except admin.sites.AlreadyRegistered:
            pass

def unregister_regulation():
    """Unregister Regulation table. Don't show any field. Gene field in the
    Regulation model ruins everything. As default behavior, Django tries to
    generate a dropdown box with __ALL__ genes in the database, which kills the
    server."""
    admin.site.unregister(Regulation)

# Register models
admin.site.register(Curator, CuratorAdmin)
admin.site.register(Curation, CurationAdmin)
admin.site.register(SiteInstance, SiteInstanceAdmin)
admin.site.register(Curation_SiteInstance, Curation_SiteInstanceAdmin)
admin.site.register(NotAnnotatedSiteInstance, NotAnnotatedSiteInstanceAdmin)
admin.site.register(Publication, PublicationAdmin)
admin.site.register(ExperimentalTechnique, ExperimentalTechniqueAdmin)
admin.site.register(Gene, GeneAdmin)
admin.site.register(Genome, GenomeAdmin)
admin.site.register(TFFamily, TFFamilyAdmin)
admin.site.register(TF, TFAdmin)
admin.site.register(NCBISubmission, NCBISubmissionAdmin)
register_rest()
unregister_regulation()

