from django.conf.urls import patterns, url


urlpatterns = patterns('baseapp.views',
    url(r'^get_weblogo/$', 'get_weblogo'),
)
