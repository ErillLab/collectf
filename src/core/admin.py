from django.contrib import admin

import models


class CurationAdmin(admin.ModelAdmin):
    filter_horizontal = ('TF_instances',)
    list_display = ('curation_id',
                    'site_species',
                    'PMID',
                    'TF',
                    'curator',
                    'validated_by',
                    'NCBI_submission_ready')
    list_filter = ('curator', 'NCBI_submission_ready', 'requires_revision')
    ordering = ('-curation_id',)


class CuratorAdmin(admin.ModelAdmin):
    list_display = ('username', 'email', 'name')


class PublicationAdmin(admin.ModelAdmin):
    list_display = ('publication_id',
                    'pmid',
                    'title',
                    'reported_species',
                    'reported_TF',
                    'assigned_to')
    list_filter = ('assigned_to', 'curation_complete')
    search_fields = ('pmid', 'reported_TF', 'reported_species')


class GeneAdmin(admin.ModelAdmin):
    list_display = ('gene_accession', 'name')
    list_filter = ('genome__genome_accession',)


class GenomeAdmin(admin.ModelAdmin):
    list_display = ('genome_accession', 'organism')


class TaxonomyAdmin(admin.ModelAdmin):
    list_display = ('taxonomy_id', 'rank', 'name')


class TFAdmin(admin.ModelAdmin):
    list_display = ('name', 'family')
    list_filter = ('family',)
    ordering = ('name',)


class TFFamilyAdmin(admin.ModelAdmin):
    list_display = ('name',)
    ordering = ('name',)


class TFInstanceAdmin(admin.ModelAdmin):
    list_display = ('uniprot_accession', 'TF', 'description')
    ordering = ('uniprot_accession', )


class SiteInstanceAdmin(admin.ModelAdmin):
    list_display = ('site_id',)
    ordering = ('-site_id', )


class Curation_SiteInstanceAdmin(admin.ModelAdmin):
    list_filter = ('site_type',)
    filter_horizontal = ('experimental_techniques',)
    ordering = ('-id',)


class NotAnnotatedSiteInstanceAdmin(admin.ModelAdmin):
    filter_horizontal = ('experimental_techniques',)
    ordering = ('-id',)


class ExperimentalTechniqueAdmin(admin.ModelAdmin):
    list_display = ('name',)
    filter_horizontal = ('categories',)


class ExperimentalTechniqueCategoryAdmin(admin.ModelAdmin):
    list_display = ('name',)


class ChipInfoAdmin(admin.ModelAdmin):
    list_display = ('chip_info_id',)


class ExternalDatabaseAdmin(admin.ModelAdmin):
    list_display = ('ext_database_id',
                    'ext_database_name',
                    'ext_database_descripton')


admin.site.register(models.Curation, CurationAdmin)
admin.site.register(models.Curator, CuratorAdmin)
admin.site.register(models.Publication, PublicationAdmin)
admin.site.register(models.Gene, GeneAdmin)
admin.site.register(models.Genome, GenomeAdmin)
admin.site.register(models.Taxonomy, TaxonomyAdmin)
admin.site.register(models.TF, TFAdmin)
admin.site.register(models.TFFamily, TFFamilyAdmin)
admin.site.register(models.TFInstance, TFInstanceAdmin)
admin.site.register(models.SiteInstance, SiteInstanceAdmin)
admin.site.register(models.Curation_SiteInstance, Curation_SiteInstanceAdmin)
admin.site.register(models.NotAnnotatedSiteInstance,
                    NotAnnotatedSiteInstanceAdmin)
admin.site.register(models.ExperimentalTechnique, ExperimentalTechniqueAdmin)
admin.site.register(models.ExperimentalTechniqueCategory,
                    ExperimentalTechniqueCategoryAdmin)
admin.site.register(models.ChipInfo, ChipInfoAdmin)
admin.site.register(models.ExternalDatabase, ExternalDatabaseAdmin)
