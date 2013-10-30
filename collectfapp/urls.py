from django.conf.urls import patterns, url

urlpatterns = patterns('collectfapp.views',
    url(r'^$', 'home'),
    url(r'^success/$', 'success'),
    url(r'^ext_pubmed_submission', 'pub_external_submission'),
)

urlpatterns += patterns('collectfapp.curationview',
    url(r'^curation/$', 'curation'),
)


urlpatterns += patterns('collectfapp.editcurationview',
    url(r'^edit_curation/(?P<cid>\d+)/$', 'edit_curation'),
)

urlpatterns += patterns('collectfapp.pubview',
    url(r'^pubmed_submission/$', 'pubmed_submission'),
    url(r'^non_pubmed_submission/$', 'non_pubmed_submission'),
)


