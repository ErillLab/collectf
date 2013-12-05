from django.conf.urls import patterns, url


urlpatterns = patterns('browseapp.json_response',
    url(r'^get_genomes/$', 'get_genomes'),
    url(r'^get_TF_instances/$', 'get_TF_instances'),
)

urlpatterns += patterns('browseapp.export',
    url(r'^export_sites/$', 'export_sites'),
)

urlpatterns += patterns('browseapp.browse_TF_and_species',
    url(r'^view_report/(?P<TF_param>\w+)/(?P<TF_ids>[\w,]+)/(?P<species_param>\w+)/(?P<species_ids>[\w,]+)/$',
        'browse_TF_and_species_selected'),
    url(r'^view_report_w_non_motif/(?P<TF_param>\w+)/(?P<TF_ids>[\w,]+)/(?P<species_param>\w+)/(?P<species_ids>[\w,]+)/$',
        'browse_TF_and_species_selected_non_motif'),
)



urlpatterns += patterns('browseapp.view_curation',
    url(r'^view_curation/(?P<cid>\d+)/$', 'view_curation'),
)

urlpatterns += patterns('browseapp.view_results',
    url(r'^view_results/$', 'view_results'),
)

urlpatterns += patterns('browseapp.view_all',
    url(r'^view_all_curations/$', 'view_all_curations'),
    url(r'^view_all_publications/$', 'view_all_publications')
)

urlpatterns += patterns('browseapp.browse_front',
    url(r'^browse_front_TF/$', 'browse_TF'),
    url(r'^browse_TF_all_reports_ajax/(?P<t>\w+)/(?P<id>\d+)/$', 'browse_TF_all_reports_ajax'),
    url(r'^browse_front_tax/$', 'browse_tax'),
    url(r'^browse_tax_all_reports_ajax/(?P<id>\d+)/$', 'browse_tax_all_reports_ajax'),
    url(r'^browse_front_tech/$', 'browse_techniques'),
    url(r'^browse_techniques_all_reports_ajax/(?P<type>\w+)/(?P<id>\d+)/$', 'browse_techniques_all_reports_ajax'),
)

urlpatterns += patterns('browseapp.search',
    url(r'^search/$', 'search')
)

urlpatterns += patterns('browseapp.compare_motifs',
    url(r'^motif_similarity_measure/$', 'motif_sim_measure', name="motif_sim_measure"),
    url(r'^compare_motifs/$', 'motif_comparison_step1', name="compare_motifs"),
    url(r'^compare_motifs/2/$', 'motif_comparison_step2'),
)

                       
