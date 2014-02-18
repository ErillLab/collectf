# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'Curation_SiteInstance.TF_function'
        db.add_column(u'base_curation_siteinstance', 'TF_function',
                      self.gf('django.db.models.fields.CharField')(default='N/A', max_length=50),
                      keep_default=False)

        # Adding M2M table for field experimental_techniques on 'Curation_SiteInstance'
        m2m_table_name = db.shorten_name(u'base_curation_siteinstance_experimental_techniques')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('curation_siteinstance', models.ForeignKey(orm[u'base.curation_siteinstance'], null=False)),
            ('experimentaltechnique', models.ForeignKey(orm[u'base.experimentaltechnique'], null=False))
        ))
        db.create_unique(m2m_table_name, ['curation_siteinstance_id', 'experimentaltechnique_id'])

        # Adding M2M table for field TF_instances on 'Curation'
        m2m_table_name = db.shorten_name(u'base_curation_TF_instances')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('curation', models.ForeignKey(orm[u'base.curation'], null=False)),
            ('tfinstance', models.ForeignKey(orm[u'base.tfinstance'], null=False))
        ))
        db.create_unique(m2m_table_name, ['curation_id', 'tfinstance_id'])


    def backwards(self, orm):
        # Deleting field 'Curation_SiteInstance.TF_function'
        db.delete_column(u'base_curation_siteinstance', 'TF_function')

        # Removing M2M table for field experimental_techniques on 'Curation_SiteInstance'
        db.delete_table(db.shorten_name(u'base_curation_siteinstance_experimental_techniques'))

        # Removing M2M table for field TF_instances on 'Curation'
        db.delete_table(db.shorten_name(u'base_curation_TF_instances'))


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'base.chipinfo': {
            'Meta': {'object_name': 'ChipInfo'},
            'assay_conditions': ('django.db.models.fields.TextField', [], {}),
            'chip_info_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_notes': ('django.db.models.fields.TextField', [], {})
        },
        u'base.curation': {
            'Meta': {'object_name': 'Curation'},
            'NCBI_submission_ready': ('django.db.models.fields.BooleanField', [], {}),
            'TF': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.TF']", 'null': 'True'}),
            'TF_function': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'TF_instance': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.TFInstance']"}),
            'TF_instances': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'curations'", 'symmetrical': 'False', 'to': u"orm['base.TFInstance']"}),
            'TF_species': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'TF_type': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'chip_info': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.ChipInfo']", 'null': 'True', 'blank': 'True'}),
            'complex_notes': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'confidence': ('django.db.models.fields.BooleanField', [], {}),
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'curation_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'curator': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curator']"}),
            'experimental_process': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'experimental_techniques': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'curation'", 'symmetrical': 'False', 'to': u"orm['base.ExperimentalTechnique']"}),
            'forms_complex': ('django.db.models.fields.BooleanField', [], {}),
            'last_modified': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'notes': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publication': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Publication']"}),
            'quantitative_data_format': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'requires_revision': ('django.db.models.fields.CharField', [], {'max_length': '20', 'null': 'True', 'blank': 'True'}),
            'site_instances': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['base.SiteInstance']", 'through': u"orm['base.Curation_SiteInstance']", 'symmetrical': 'False'}),
            'site_species': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'validated_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'validated_by'", 'null': 'True', 'to': u"orm['base.Curator']"})
        },
        u'base.curation_externaldatabase': {
            'Meta': {'object_name': 'Curation_ExternalDatabase'},
            'accession_number': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'curation': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curation']"}),
            'external_database': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.ExternalDatabase']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'base.curation_siteinstance': {
            'Meta': {'object_name': 'Curation_SiteInstance'},
            'TF_function': ('django.db.models.fields.CharField', [], {'default': "'N/A'", 'max_length': '50'}),
            'annotated_seq': ('django.db.models.fields.TextField', [], {'max_length': '100000'}),
            'curation': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curation']"}),
            'experimental_techniques': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['base.ExperimentalTechnique']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_obsolete': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'quantitative_value': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'regulates': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['base.Gene']", 'through': u"orm['base.Regulation']", 'symmetrical': 'False'}),
            'site_instance': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.SiteInstance']"}),
            'site_type': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'why_obsolete': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'})
        },
        u'base.curator': {
            'Meta': {'object_name': 'Curator'},
            'curator_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'curator_type': ('django.db.models.fields.CharField', [], {'default': "'external'", 'max_length': '20'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        },
        u'base.experimentaltechnique': {
            'Meta': {'object_name': 'ExperimentalTechnique'},
            'categories': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['base.ExperimentalTechniqueCategory']", 'symmetrical': 'False'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'preset_function': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True'}),
            'technique_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'base.experimentaltechniquecategory': {
            'Meta': {'object_name': 'ExperimentalTechniqueCategory'},
            'category_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'base.externaldatabase': {
            'Meta': {'object_name': 'ExternalDatabase'},
            'ext_database_descripton': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'ext_database_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ext_database_name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '50'}),
            'ext_database_url_format': ('django.db.models.fields.CharField', [], {'max_length': '500'})
        },
        u'base.gene': {
            'Meta': {'object_name': 'Gene'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '1000'}),
            'end': ('django.db.models.fields.IntegerField', [], {}),
            'gene_accession': ('django.db.models.fields.CharField', [], {'max_length': '30', 'null': 'True', 'blank': 'True'}),
            'gene_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Genome']"}),
            'locus_tag': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'start': ('django.db.models.fields.IntegerField', [], {}),
            'strand': ('django.db.models.fields.IntegerField', [], {})
        },
        u'base.genome': {
            'GC_content': ('django.db.models.fields.FloatField', [], {}),
            'Meta': {'object_name': 'Genome'},
            'chromosome': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'genome_accession': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'genome_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'genome_sequence': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['base.GenomeSequence']", 'unique': 'True'}),
            'gi': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'organism': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'taxonomy': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Taxonomy']"})
        },
        u'base.genomesequence': {
            'Meta': {'object_name': 'GenomeSequence'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {})
        },
        u'base.ncbisubmission': {
            'Meta': {'object_name': 'NCBISubmission'},
            'curation_site_instance': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curation_SiteInstance']"}),
            'genome_submitted_to': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'submission_time': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'})
        },
        u'base.notannotatedsiteinstance': {
            'Meta': {'object_name': 'NotAnnotatedSiteInstance'},
            'curation': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curation']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'max_length': '100000'})
        },
        u'base.publication': {
            'Meta': {'object_name': 'Publication'},
            'assigned_to': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curator']", 'null': 'True', 'blank': 'True'}),
            'authors': ('django.db.models.fields.CharField', [], {'max_length': '1000'}),
            'contains_expression_data': ('django.db.models.fields.BooleanField', [], {}),
            'contains_promoter_data': ('django.db.models.fields.BooleanField', [], {}),
            'curation_complete': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'issue': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'journal': ('django.db.models.fields.CharField', [], {'max_length': '1000'}),
            'pages': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'pdf': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'pmid': ('django.db.models.fields.CharField', [], {'max_length': '30', 'null': 'True', 'blank': 'True'}),
            'publication_date': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'publication_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'publication_type': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'reported_TF': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'reported_species': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'submission_notes': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '1000'}),
            'url': ('django.db.models.fields.CharField', [], {'max_length': '1000', 'null': 'True', 'blank': 'True'}),
            'volume': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'base.regulation': {
            'Meta': {'object_name': 'Regulation'},
            'curation_site_instance': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Curation_SiteInstance']"}),
            'evidence_type': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'base.siteinstance': {
            'Meta': {'object_name': 'SiteInstance'},
            '_seq': ('django.db.models.fields.TextField', [], {'max_length': '100000'}),
            'end': ('django.db.models.fields.IntegerField', [], {}),
            'genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Genome']"}),
            'site_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'start': ('django.db.models.fields.IntegerField', [], {}),
            'strand': ('django.db.models.fields.IntegerField', [], {})
        },
        u'base.taxonomy': {
            'Meta': {'object_name': 'Taxonomy'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.Taxonomy']", 'null': 'True'}),
            'rank': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'taxonomy_id': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'})
        },
        u'base.tf': {
            'Meta': {'object_name': 'TF'},
            'TF_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'family': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['base.TFFamily']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'base.tffamily': {
            'Meta': {'object_name': 'TFFamily'},
            'TF_family_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'base.tfinstance': {
            'Meta': {'object_name': 'TFInstance'},
            'description': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'protein_accession': ('django.db.models.fields.CharField', [], {'max_length': '20', 'primary_key': 'True'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['base']