from collectfapp.models import *
from django.contrib import admin
from django.db.models import get_models
from django.db.models import get_app

# curation admin view
class CurationAdmin(admin.ModelAdmin):
    def PMID(self, obj):
        return obj.publication.pmid
    
    filter_horizontal = ("experimental_techniques",)
    list_display = ('curation_id', 'TF_species', 'PMID', 'curator')
    list_filter = ('curator', 'NCBI_submission_ready', 'requires_revision')
    ordering = ('-curation_id',)

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
    list_filter = ('genome__genome_accession',)
    ordering = ('-site_id', )
    
class Curation_SiteInstanceAdmin(admin.ModelAdmin):
    list_display = ('id',)
    raw_id_fields = ('chipseq_info',)
    list_filter = ('is_motif_associated',)
    ordering = ('-id',)
    
class PublicationAdmin(admin.ModelAdmin):
    list_display = ('publication_id', 'pmid', 'title', 'assigned_to')
    list_filter = ('assigned_to', 'curation_complete', 'reported_TF', 'reported_species')

class ExperimentalTechniqueAdmin(admin.ModelAdmin):
    list_display = ('name', )

class GeneAdmin(admin.ModelAdmin):
    list_display = ('gene_accession', 'name')
    list_filter = ('genome__genome_accession',)

class GenomeAdmin(admin.ModelAdmin):
    def strain_name(self, obj):
        return obj.strain.name
    list_display = ('genome_accession', 'strain_name')

class TFAdmin(admin.ModelAdmin):
    list_display = ('TF_id', 'name', 'family')
    list_filter = ('family',)

# register edited models
admin.site.register(Curation, CurationAdmin)
admin.site.register(Curator, CuratorAdmin)
admin.site.register(SiteInstance, SiteInstanceAdmin)
admin.site.register(Curation_SiteInstance, Curation_SiteInstanceAdmin)
admin.site.register(Publication, PublicationAdmin)
admin.site.register(ExperimentalTechnique, ExperimentalTechniqueAdmin)
admin.site.register(Gene, GeneAdmin)
admin.site.register(Genome, GenomeAdmin)
admin.site.register(TF, TFAdmin)


# register rest of models
app = get_app("collectfapp")
for model in get_models(app):
    try:
        admin.site.register(model)
    except admin.sites.AlreadyRegistered:
        pass

# unregister Regulation
# don't show any field. Gene field in the Regulation model
# ruins everything. As default behavior, Django tries to
# generate a dropdown box with __ALL__ genes in the database,
# which kills the server.
admin.site.unregister(Regulation)
