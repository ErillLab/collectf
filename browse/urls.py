from django.conf.urls import patterns, url

urlpatterns = patterns(
    'browse.browse_tax',
    (r'browse_tax/$', 'browse_taxonomy'),
    (r'get_results_tax/(?P<object_id>\d+)$', 'get_results_taxonomy'),
)

urlpatterns += patterns(
    'browse.browse_TF',
    (r'browse_tf/$', 'browse_TF'),
    (r'get_results_TF/(?P<object_id>\d+)/$', 'get_results_TF'),
    (r'get_results_TF_family/(?P<object_id>\d+)/$', 'get_results_TF_family'),
)

urlpatterns += patterns(
    'browse.browse_tech',
    (r'browse_tech/$', 'browse_tech'),
    (r'get_results_technique/(\d+)$', 'get_results_technique'),
    (r'get_results_technique_all/(binding|expression)$', 'get_results_all'),
    (r'get_results_technique_category/(binding|expression)/(\d+)$',
     'get_results_category'),
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
    (r'view_reports_by_TF_and_species/(\d+)/(\d+)/$',
     'view_reports_by_TF_and_species'),
    (r'view_reports_by_id_list', 'view_reports_by_id_list'),
    (r'view_reports_by_TF_family/(\d+)$', 'view_reports_by_TF_family'),
    (r'view_reports_by_TF/(\d+)$', 'view_reports_by_TF'),
    (r'view_reports_by_taxonomy/(\d+)$', 'view_reports_by_taxonomy'),
    (r'view_reports_by_all_techniques/(binding|expression)$',
     'view_reports_by_all_techniques'),
    (r'view_reports_by_technique_category/(binding|expression)/(\d+)$',
     'view_reports_by_technique_category'),
    (r'view_reports_by_technique/(\d+)$', 'view_reports_by_technique'),
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
    url(r'^tf_list/$', 'list_tfs'),
    url(r'^species_list/$', 'list_species'),
    url(r'^exp_tech_list/$', 'list_experimental_techniques'),
    url(r'^release_history/$', 'release_history'),
    url(r'^view_all_curations/$', 'view_all_curations'),
    url(r'^view_all_publications/$', 'view_all_publications'),
    url(r'^update_stats/$', 'update_stats'),
    url(r'^list_all_motifs/$', 'list_all_motifs'),
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
