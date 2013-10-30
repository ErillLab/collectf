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

    url(r'^admin/', include(admin.site.urls)),
    url(r'^collectfapp/', include('collectfapp.urls')),
    url(r'^browseapp/', include('browseapp.urls')),
    url(r'^baseapp/', include('baseapp.urls')),
    url(r'^mainpageapp/', include('mainpageapp.urls')),
    url(r'^ncbiapp/', include('ncbiapp.urls')),
    url(r'^dbstatsapp/', include('dbstatsapp.urls')),
    url(r'^mastercuratorapp/', include('mastercuratorapp.urls')),
)
urlpatterns += staticfiles_urlpatterns()
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

