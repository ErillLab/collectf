from django.conf.urls import patterns, include, url
from collectfapp.signupview import *

from django.views.generic.simple import direct_to_template


import collectfapp.views
import collectfapp.signupview as signupview
import collectfapp.pubview as pubview
import collectfapp.curationview as curationview
import collectfapp.editcurationview as editcurationview

import browseapp.browse_TF_and_species
import browseapp.search
import browseapp.compare_motifs
import browseapp.view_results
import browseapp.view_all
import browseapp.view_site
import browseapp.view_curation
import browseapp.json_response
import browseapp.browse_front
import browseapp.export
import dbstatsapp.views
import ncbiapp.views
import baseapp.views

import mainpageapp.views

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
    url(r'^accounts/register/$', mainpageapp.views.register_request),
    # login
    url(r'^accounts/login/$', signupview.login),
    # logout
    url(r'^accounts/logout/$', signupview.logout),
    # main page
    url(r'^$', collectfapp.views.home),

    # external pubmed submission
    url(r'^ext_pubmed_submission', collectfapp.views.pub_external_submission),
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
    url(r'^view_all_curations/$', browseapp.view_all.view_all_curations),
    # view all publications
    url(r'^view_all_publications/$', browseapp.view_all.view_all_publications),
    # search
    url(r'^search/$', browseapp.search.search),

    # motif comparison
    url(r'^compare_motifs/$', browseapp.compare_motifs.MotifComparisonWizard.as_view(browseapp.compare_motifs.FORMS), name="compare_motifs"),
                       
    # view results
    url(r'^view_results/$', browseapp.view_results.view_results),
     
    # browse by taxonomy
    url(r'^browse_front_TF/$', browseapp.browse_front.browse_TF),
    url(r'^browse_TF_all_reports_ajax/(?P<t>\w+)/(?P<id>\d+)/$', browseapp.browse_front.browse_TF_all_reports_ajax),
    url(r'^browse_front_tax/$', browseapp.browse_front.browse_tax),                   
    url(r'^browse_tax_all_reports_ajax/(?P<id>\d+)/$', browseapp.browse_front.browse_tax_all_reports_ajax),
    url(r'^browse_front_tech/$', browseapp.browse_front.browse_techniques),
    url(r'^browse_techniques_all_reports_ajax/(?P<type>\w+)/(?P<id>\d+)/$', browseapp.browse_front.browse_techniques_all_reports_ajax),
    # view curation
    url(r'^view_curation/(?P<cid>\d+)/$', browseapp.view_curation.view_curation),
    # browse site
    url(r'^expsite_(?P<dbxref_id>\w+)$', browseapp.view_site.browse_by_site),
    url(r'^EXPSITE_(?P<dbxref_id>\w+)$', browseapp.view_site.browse_by_site),
    # view report
    url(r'^view_report/(?P<TF_param>\w+)/(?P<TF_ids>[\w,]+)/(?P<species_param>\w+)/(?P<species_ids>[\w,]+)/$',
        browseapp.browse_TF_and_species.browse_TF_and_species_selected),
    # view report (integrate non-motif-associated sites)
    url(r'^view_report_w_non_motif/(?P<TF_param>\w+)/(?P<TF_ids>[\w,]+)/(?P<species_param>\w+)/(?P<species_ids>[\w,]+)/$',
        browseapp.browse_TF_and_species.browse_TF_and_species_selected_non_motif),
           
    # export fasta/csv
    url(r'^export_sites/$', browseapp.export.export_sites),
    # database statistics
    url(r'^db_stats/$', dbstatsapp.views.curation_stats),
    url(r'^curator_roster/$', dbstatsapp.views.curator_roster),
    url(r'^exp_tech_list/$', dbstatsapp.views.experimental_techniques),
    # export tbl for ncbi submission
    url(r'^export_ncbi/$', ncbiapp.views.export_tbl_view),

    # JSON requests
    url(r'^get_genomes/$', browseapp.json_response.get_genomes),
    url(r'^get_TF_instances/$', browseapp.json_response.get_TF_instances),
    url(r'^get_weblogo', baseapp.views.get_weblogo),

                       
    # main page handlers
    url(r'^main_page_greet/$', mainpageapp.views.greet),
    url(r'^main_page_about/$', mainpageapp.views.about),
    url(r'^main_page_browse/$', mainpageapp.views.browse),
    url(r'^main_page_search/$', mainpageapp.views.search),
    url(r'^main_page_compare_motifs/$', mainpageapp.views.compare_motifs),
    url(r'^main_page_contribute/$', mainpageapp.views.contribute),
    url(r'^main_page_feedback/$', mainpageapp.views.feedback),
    url(r'^main_page_stats/$', mainpageapp.views.stats),
    url(r'^main_page_cite/$', mainpageapp.views.cite),
    url(r'^main_page_links/$', mainpageapp.views.links),
    url(r'^main_page_ack/$', mainpageapp.views.acknowledgements),
    url(r'^feedback_send_email/$', mainpageapp.views.feedback_send_email),
                       
                       
)

from django.contrib.staticfiles.urls import staticfiles_urlpatterns
urlpatterns += staticfiles_urlpatterns()
