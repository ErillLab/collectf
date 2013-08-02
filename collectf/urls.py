from django.conf.urls import patterns, include, url
from collectfapp.signupview import *

from django.views.generic.simple import direct_to_template


import collectfapp.views
import collectfapp.signupview as signupview
import collectfapp.pubview as pubview
import collectfapp.curationview as curationview
import collectfapp.editcurationview as editcurationview

import browseapp.browse_species
import browseapp.browse_TF
import browseapp.browse_TF_and_species
import browseapp.browse_all
import browseapp.browse_site
import browseapp.browse_curation
import browseapp.json_response
import dbstatsapp.views
import ncbiapp.views

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'collectf.views.home', name='home'),
    # url(r'^collectf/', include('collectf.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
    # registration
    # registration is not allowed, admin creates new accounts.
    # url(r'^accounts/register/$', signupview.register),
    # login
    url(r'^accounts/login/$', signupview.login),
    # logout
    url(r'^accounts/logout/$', signupview.logout),
    # main page
    url(r'^$', collectfapp.views.home),
    # pubmed publication submission
    url(r'^pubmed_submission/$', pubview.pubmed_submission),
    # nonpubmed publication submission
    url(r'^non_pubmed_submission/$', pubview.non_pubmed_submission),
    # curation
    url(r'^curation/$', curationview.curation),
    # edit curation
    url(r'^edit_curation/(?P<cid>\d+)/$', editcurationview.edit_curation),
                       
    # success page
    url(r'^success/$', collectfapp.views.success),
    # view all curations
    url(r'^view_all_curations/$', browseapp.browse_all.view_all_curations),
    # view all publications
    url(r'^view_all_publications/$', browseapp.browse_all.view_all_publications),
    # browse
    url(r'^browse/$', browseapp.browse_TF_and_species.browse_TF_and_species),
    # browse by TF and species
    url(r'^browse_TF_sp/(?P<TF_id>\d+)/(?P<species_id>\d+)/$', browseapp.browse_TF_and_species.browse_TF_and_species_selected),
    # browse by TF main
    url(r'^browse_tf_main/$', browseapp.browse_TF.browse_by_TF_main),
    # browse by TF family
    url(r'^browse_tf_family/(?P<TF_family_id>\d+)$', browseapp.browse_TF.browse_by_TF_family),
    # browse by TF
    url(r'^browse_tf/(?P<TF_id>\d+)$', browseapp.browse_TF.browse_by_TF),
    # browse by species main
    url(r'^browse_sp_main/$', browseapp.browse_species.browse_by_species_main),
    # browse by taxon elms
    url(r'^browse_sp_taxon/(?P<tax_id>\d+)/$', browseapp.browse_species.browse_by_species_taxon),
    # browse by species
    url(r'^browse_sp/(?P<sp_tax_id>\d+)/$', browseapp.browse_species.browse_by_species),

                       
    # browse curation
    url(r'^browse_curation/(?P<cid>\d+)/$', browseapp.browse_curation.browse_curation),
    # browse site
    url(r'^expsite_(?P<dbxref_id>\w+)$', browseapp.browse_site.browse_by_site),
    url(r'^EXPSITE_(?P<dbxref_id>\w+)$', browseapp.browse_site.browse_by_site),
           
    # export fasta/csv
    url(r'^export_sites/$', browseapp.browse_TF_and_species.export_sites),
    # database statistics
    url(r'^db_stats/$', dbstatsapp.views.curation_stats),
    # export tbl for ncbi submission
    url(r'^export_ncbi/$', ncbiapp.views.export_tbl_view),

    # JSON requests
    url(r'^get_genomes/$', browseapp.json_response.get_genomes),
    url(r'^get_TF_instances/$', browseapp.json_response.get_TF_instances),

                       
)

from django.contrib.staticfiles.urls import staticfiles_urlpatterns
urlpatterns += staticfiles_urlpatterns()
