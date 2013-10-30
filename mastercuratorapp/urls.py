from django.conf.urls import patterns, url


urlpatterns = patterns('mastercuratorapp.views',
    url(r'^home/$', 'home'),
    url(r'^validate_curation/(?P<curation_id>\d+)$', 'validate_curation'),
    url(r'^edit_curation/(?P<curation_id>\d+)$', 'edit_curation'),
)
