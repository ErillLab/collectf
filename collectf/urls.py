from django.conf.urls import patterns, include, url
from collectfapp.signup import *
from collectfapp.views import *
from django.contrib.auth.decorators import login_required

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
    url(r'^accounts/register/$', register),
    # login
    url(r'^accounts/login/$', login, {'extra_context': {'next': '/'}}),
    # logout
    url(r'^accounts/logout/$', logout),
    # main page
    url(r'^$', login_required(home_view)),
    
)

print login
