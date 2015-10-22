
from django.conf.urls import patterns, include, url

from django.contrib import admin
admin.autodiscover()

import base.views
import browse.view_site
import browse.view_reports

urlpatterns = patterns('',
    url(r'^admin/', include(admin.site.urls)),
    url(r'', include('base.urls')),
    url(r'^home/', include('homepage.urls')),
    url(r'^browse/', include('browse.urls')),
    url(r'^curate/', include('curate.urls')),
    url(r'^ncbi/', include('ncbi.urls')),

    # login/logout
    url(r'^accounts/login/$', base.views.login),
    url(r'^accounts/logout/$', base.views.logout),
    url(r'^accounts/', include('registration.backends.default.urls')),

    # This should ideally be in browseapp urls.py, but since NCBI links cannot
    # be changed, it will serve here.
    url(r'^expsite_(?P<dbxref_id>\w+)$', browse.view_site.view_site),
    url(r'^EXPSITE_(?P<dbxref_id>\w+)$', browse.view_site.view_site),
    # UniProt link-outs
    url(r'^expreg_(?P<uniprot_dbxref>\w+)$',
        browse.view_reports.view_reports_by_uniprot_dbxref),
    url(r'^EXPREG_(?P<uniprot_dbxref>\w+)$',
        browse.view_reports.view_reports_by_uniprot_dbxref),
    url(r'^uniprot_(?P<uniprot_accession>\w+)$',
        browse.view_reports.view_reports_by_uniprot_accession),
)
