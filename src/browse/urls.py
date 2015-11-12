from django.conf.urls import url

from . import test_view
from . import homepage_views
from . import browse_by_TF
from . import browse_by_taxonomy
from . import browse_by_technique
from . import view_motif_reports
from . import search_motifs
from . import compare_motifs
from . import view_curation
from . import list_all

urlpatterns = [
    url(r'^test_view$', test_view.test_view),
    # Homepage view
    url(r'^home/$', homepage_views.home, name='homepage_home'),

    url(r'^about/$', homepage_views.about, name='homepage_about'),

    url(r'^browse/$', homepage_views.browse, name='homepage_browse'),

    url(r'^search/$', homepage_views.search, name='homepage_search'),

    url(r'^compare/$', homepage_views.compare, name='homepage_compare'),

    url(r'^contribute/$', homepage_views.contribute,
        name='homepage_contribute'),

    url(r'^stats/$', homepage_views.stats, name='homepage_stats'),

    url(r'^links/$', homepage_views.links, name='homepage_links'),

    url(r'^feedback/$', homepage_views.feedback, name='homepage_feedback'),

    url(r'^cite/$', homepage_views.cite, name='homepage_cite'),

    url(r'^acknowledgments/$', homepage_views.acknowledgments,
        name='homepage_acknowledgments'),

    # Browse by TF
    url(r'^browse_by_TF/$', browse_by_TF.browse_TF, name='browse_by_TF'),

    url(r'^get_results_by_TF_family/(?P<object_id>\d+)/$',
        browse_by_TF.get_results_by_TF_family),

    url(r'^get_results_by_TF/(?P<object_id>\d+)/$',
        browse_by_TF.get_results_by_TF),

    # Browse by taxonomy
    url(r'^browse_by_taxonomy/$', browse_by_taxonomy.browse_taxonomy,
        name='browse_by_taxonomy'),

    url(r'^get_results_by_taxonomy/(?P<object_id>\d+)/$',
        browse_by_taxonomy.get_results_by_taxonomy),

    # Browse by experimental technique
    url(r'^browse_by_technique/$', browse_by_technique.browse_technique,
        name='browse_by_technique'),

    url(r'get_results_by_technique/(\d+)/$',
        browse_by_technique.get_results_by_technique),

    url(r'get_results_by_technique_function/(binding|expression)/$',
        browse_by_technique.get_results_by_technique_function),

    url(r'get_results_by_technique_category/(binding|expression)/(\d+)/$',
        browse_by_technique.get_results_by_technique_category),

    # View motif reports by TF and species
    url(r'^view_motif_reports_by_TF_and_species/(\d+)/(\d+)/$',
        view_motif_reports.view_reports_by_TF_and_species,
        name='view_motif_reports_by_TF_and_species'),

    # View motif reports by TF family
    url(r'^view_motif_reports_by_TF_family/(\d+)/$',
        view_motif_reports.view_reports_by_TF_family),

    # View motif reports by TF
    url(r'^view_motif_reports_by_TF/(\d+)/$',
        view_motif_reports.view_reports_by_TF),

    # View motif reports by technique function
    url(r'^view_motif_reports_by_technique_function/(binding|expression)/$',
        view_motif_reports.view_reports_by_technique_function),

    # View motif reports by technique category
    url(r'^view_motif_reports_by_technique_category/(binding|expression)/(\d+)/$',
        view_motif_reports.view_reports_by_technique_category),

    # View motif reports by technique
    url(r'^view_motif_reports_by_technique/(\d+)/$',
        view_motif_reports.view_reports_by_technique),

    # View motif reports by taxonomy
    url(r'^view_motif_reports_by_taxonomy/(\d+)/$',
        view_motif_reports.view_reports_by_taxonomy),

    # View customized motif reports by list of Curation_SiteInstance IDs.
    url(r'^view_motif_reports_by_curation_site_instance_ids/$',
        view_motif_reports.view_reports_by_curation_site_instance_ids,
        name='view_reports_by_id'),

    # Search
    url(r'^search_motifs/$', search_motifs.search, name='motif_search'),

    # Motif comparison
    url(r'^compare_motifs/$',
        compare_motifs.motif_comparison_search_first_motif,
        name='motif_comparison_first_step'),
    url(r'^compare_motifs/2/$',
        compare_motifs.motif_comparison_search_second_motif,
        name='motif_comparison_second_step'),
    url(r'^motif_similarity_measure/$',
        compare_motifs.motif_similarity_measure,
        name='motif_similarity_measure'),

    # View curation
    url(r'^view_curation/(\d+)/$', view_curation.view_curation,
        name='view_curation'),

    # List all TFs
    url(r'^list_all_TFs/$', list_all.list_all_TFs,
        name='list_all_TFs'),

    # List all species
    url(r'^list_all_species/$', list_all.list_all_species,
        name='list_all_species'),

    # List all experimental techniques
    url(r'^list_all_experimental_techniques/$',
        list_all.list_all_experimental_techniques,
        name='list_all_experimental_techniques'),

    # List all publications
    url(r'^list_all_publications/$', list_all.list_all_publications,
        name='list_all_publications'),

    # List all curations
    url(r'^list_all_curations/$', list_all.list_all_curations,
        name='list_all_curations'),
    

]
