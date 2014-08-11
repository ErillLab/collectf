from django.conf.urls import patterns, url

urlpatterns = patterns(
    'browse.browse_tax',
    (r'browse_tax/$', 'browse_taxonomy'),
    (r'get_results_tax/(?P<taxid>\d+)$', 'get_results_taxonomy'),
)

urlpatterns += patterns(
    'browse.browse_TF',
    (r'browse_tf/$', 'browse_TF'),
    (r'get_results_tf/(?P<type_>\w+)/(?P<id_>\d+)/$', 'get_results_TF'),
)

urlpatterns += patterns(
    'browse.browse_tech',
    (r'browse_tech/$', 'browse_tech'),
    url(r'^get_results_tech/(?P<type_>\w+)/(?P<id_>\d+)/$', 'get_results_tech'),
)

urlpatterns += patterns(
    'browse.search',
    (r'search/$', 'search'),
)

urlpatterns += patterns(
    'browse.view_curation',
    (r'^view_curation/(?P<cid>\d+)/$', 'view_curation'),
)

urlpatterns += patterns(
    'browse.view_reports',
    (r'view_reports/(?P<tax_param_type>.+)/(?P<tax_param>-?\d+)/(?P<tf_param_type>.+)/(?P<tf_param>-?\d+)/(?P<tech_param_type>.+)/(?P<tech_param>-?\d+)/(?P<integrate_non_motif>\d)',
     'view_reports'),
)

urlpatterns += patterns(
    'browse.json_views',
    (r'get_genomes/$', 'get_genomes'),
    (r'get_TF_instances/$', 'get_TF_instances'),
)

urlpatterns += patterns(
    'browse.stats',
    url(r'^curation_stats/$', 'curation_stats'),
    url(r'^curator_roster/$', 'curator_roster'),
    url(r'^exp_tech_list/$', 'experimental_techniques'),
    url(r'^release_history/$', 'release_history'),
    url(r'^view_all_curations/$', 'view_all_curations'),
    url(r'^view_all_publications/$', 'view_all_publications'),
    url(r'^update_stats/$', 'update_stats'),
)

urlpatterns += patterns(
    'browse.export',
    url(r'^export/$', 'export_sites'),
)

urlpatterns += patterns(
    'browse.compare',
    url(r'^motif_similarity_measure/$', 'motif_sim_measure', name="motif_sim_measure"),
    url(r'^compare_motifs/$', 'motif_comparison_step1', name="compare_motifs"),
    url(r'^compare_motifs/2/$', 'motif_comparison_step2'),
)
