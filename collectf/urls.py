from django.conf.urls import patterns, include, url
from collectfapp.signupview import *

from django.views.generic.simple import direct_to_template


import collectfapp.views
import collectfapp.signupview as signupview
import collectfapp.pubview as pubview
import collectfapp.curationview as curationview
import collectfapp.editcurationview as editcurationview
import collectfapp.viewcurationview as viewcurationview
import browseapp.views

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
    # view curation
    url(r'^view_curation/(?P<cid>\d+)/$', viewcurationview.view_curation),
    # success page
    url(r'^success/$', collectfapp.views.success),
    # view all curations
    url(r'^view_all_curations/$', browseapp.views.view_all_curations),
    # view all publications
    url(r'^view_all_publications/$', browseapp.views.view_all_publications),
    # browse
    url(r'^browse/$', browseapp.views.browse),
    # browse by TF and species
    url(r'^browse_TF_sp/(?P<TF_id>\d+)/(?P<species_id>\d+)/$', browseapp.views.browse_post_TF_sp),
    # browse by species
    url(r'^browse_sp_main/$', browseapp.views.browse_by_species_main),
    # browse by TF main
    url(r'^browse_tf_main/$', browseapp.views.browse_by_TF_main),
    # browse by TF family
    url(r'^browse_tf_family/(?P<TF_family_id>\d+)$', browseapp.views.browse_by_TF_family),
    # browse by TF
    url(r'^browse_tf/(?P<TF_id>\d+)$', browseapp.views.browse_by_TF),
    # download sites (FASTA format)
    url(r'^motif_download/(?P<TF_id>\d+)/(?P<species_id>\d+)/$', browseapp.views.report_FASTA),
    # database statistics
    url(r'^db_stats/$', browseapp.views.curation_stats),
    # display weblogo
    url(r'^display_weblogo/$', browseapp.views.display_weblogo),
            
)

from django.contrib.staticfiles.urls import staticfiles_urlpatterns
urlpatterns += staticfiles_urlpatterns()
