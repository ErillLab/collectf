from django.conf.urls import patterns, include, url
from django.contrib import admin
admin.autodiscover()
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf.urls.static import static
from django.conf import settings

import collectfapp.views
import collectfapp.signupview

urlpatterns = patterns('',
    url(r'^$', collectfapp.views.home),
    url(r'^accounts/login/$', collectfapp.signupview.login),
    url(r'^accounts/logout/$', collectfapp.signupview.logout),
    url(r'^accounts/', include('registration.backends.default.urls')),

    
    url(r'^admin/', include(admin.site.urls)),
    url(r'^collectfapp/', include('collectfapp.urls')),
    url(r'^browseapp/', include('browseapp.urls')),
    url(r'^baseapp/', include('baseapp.urls')),
    url(r'^mainpageapp/', include('mainpageapp.urls')),
    url(r'^ncbiapp/', include('ncbiapp.urls')),
    url(r'^dbstatsapp/', include('dbstatsapp.urls')),
    url(r'^mastercuratorapp/', include('mastercuratorapp.urls')),
)

import browseapp.view_site
# This should be in browseapp urls.py, but since NCBI links cannot be changed, it
# will serve here.
urlpatterns += patterns('',
    url(r'^expsite_(?P<dbxref_id>\w+)$', browseapp.view_site.browse_by_site),
    url(r'^EXPSITE_(?P<dbxref_id>\w+)$', browseapp.view_site.browse_by_site),
)
urlpatterns += staticfiles_urlpatterns()
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

