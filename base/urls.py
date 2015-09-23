from django.conf.urls import patterns, url

urlpatterns = patterns('base.views',
    url(r'^$', 'home'),
    url(r'^get_weblogo/$', 'get_weblogo'),
    url(r'^login$', 'login'),
    url(r'^logout$', 'logout'),
)
