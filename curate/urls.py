from django.conf.urls import patterns, url

urlpatterns = patterns(
    'curate.add_publication',
    url(r'^add_pubmed/$', 'pubmed_submission'),
    url(r'^add_non_pubmed/$', 'non_pubmed_submission'),
)

urlpatterns += patterns(
    'curate.add_curation',
    url(r'^curation/$', 'curation'),
    url(r'^high_throughput_curation/$', 'high_throughput_curation'),
)

urlpatterns += patterns(
    'curate.edit_curation',
    url(r'^edit_curation/(?P<cid>\d+)/$', 'edit_curation'),
)

urlpatterns += patterns(
    'curate.add_TF',
    url(r'^add_TF/$', 'add_TF'),
    url(r'^add_TF_family/$', 'add_TF_family'),
)

urlpatterns += patterns(
    'curate.add_technique',
    url(r'^add_technique/$', 'add_technique'),
)

urlpatterns += patterns(
    'curate.validate_curation',
    url(r'^validate_home/$', 'home'),
    url(r'^validate_curation/(?P<curation_id>\d+)$', 'validate_curation'),
)
