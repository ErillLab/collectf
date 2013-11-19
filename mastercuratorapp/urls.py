from django.conf.urls import patterns, url

urlpatterns = patterns('mastercuratorapp.views',
    url(r'^home/$', 'home'),
    url(r'^validate_edit_main/$', 'validate_edit_main'),
    url(r'^validate_curation/(?P<curation_id>\d+)$', 'validate_curation'),
    url(r'^edit_curation/(?P<curation_id>\d+)$', 'edit_curation'),
    url(r'^view_validated_curations/$', 'view_validated_curations'),
    url(r'^edit_validated_curation/(?P<curation_id>\d+)$', 'edit_validated_curation'),

)
