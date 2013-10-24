from django.conf.urls import patterns, url


urlpatterns = patterns('dbstatsapp.views',
    url(r'^curation_stats/$', 'curation_stats'),
    url(r'^curator_roster/$', 'curator_roster'),
    url(r'^exp_tech_list/$', 'experimental_techniques'),
    url(r'^release_history/$', 'release_history'),
)
