from django.conf.urls import patterns, include, url
from collectfapp.signupview import *
from django.contrib.auth.decorators import login_required

import collectfapp.views as views
import collectfapp.signupview as signupview
import collectfapp.pubview as pubview
import collectfapp.curationview as curationview
import collectfapp.editcurationview as editcurationview

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
    url(r'^accounts/register/$', signupview.register),
    # login
    url(r'^accounts/login/$', signupview.login, {'extra_context': {'next': '/'}}),
    # logout
    url(r'^accounts/logout/$', signupview.logout),
    # main page
    url(r'^$', login_required(views.home)),
    # pubmed publication submission
    url(r'^pubmed_submission/$', login_required(pubview.pubmed_submission)),
    # nonpubmed publication submission
    url(r'^non_pubmed_submission/$', login_required(pubview.non_pubmed_submission)),
    # curation
    url(r'^curation/$', login_required(curationview.curation)),
    # edit curation
    url(r'^edit_curation/(?P<cid>\d+)/$', login_required(editcurationview.edit_curation)),
)

print login
