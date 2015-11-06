from django.conf.urls import url

from . import test_view
from . import homepage_views
from . import browse_by_TF
from . import browse_by_taxonomy
from . import browse_by_technique
from . import view_motif_reports

urlpatterns = [
    url(r'^test_view$', test_view.test_view),
    # Homepage view
    url(r'^home$', homepage_views.home, name='homepage_home'),

    url(r'^about$', homepage_views.about, name='homepage_about'),

    url(r'^browse$', homepage_views.browse, name='homepage_browse'),

    url(r'^search$', homepage_views.search, name='homepage_search'),

    url(r'^compare$', homepage_views.compare, name='homepage_compare'),

    url(r'^contribute$', homepage_views.contribute,
        name='homepage_contribute'),

    url(r'^stats$', homepage_views.stats, name='homepage_stats'),

    url(r'^links$', homepage_views.links, name='homepage_links'),

    url(r'^feedback$', homepage_views.feedback, name='homepage_feedback'),

    url(r'^cite$', homepage_views.cite, name='homepage_cite'),

    url(r'^acknowledgements$', homepage_views.acknowledgements,
        name='homepage_acknowledgements'),

    # Browse by TF
    url(r'^browse_by_TF$', browse_by_TF.browse_TF, name='browse_by_TF'),

    url(r'^get_results_by_TF_family/(?P<object_id>\d+)$',
        browse_by_TF.get_results_by_TF_family),

    url(r'^get_results_by_TF/(?P<object_id>\d+)$',
        browse_by_TF.get_results_by_TF),

    # Browse by taxonomy
    url(r'^browse_by_taxonomy$', browse_by_taxonomy.browse_taxonomy,
        name='browse_by_taxonomy'),

    url(r'^get_results_by_taxonomy/(?P<object_id>\d+)$',
        browse_by_taxonomy.get_results_by_taxonomy),

    # Browse by experimental technique
    url(r'^browse_by_technique$', browse_by_technique.browse_technique,
        name='browse_by_technique'),

    url(r'get_results_by_technique/(\d+)/$',
        browse_by_technique.get_results_by_technique),

    url(r'get_results_by_technique_function/(binding|expression)/$',
        browse_by_technique.get_results_by_technique_function),

    url(r'get_results_by_technique_category/(binding|expression)/(\d+)/$',
        browse_by_technique.get_results_by_technique_category),

    # View motif reports
    url(r'^view_motif_reports_by_TF_family/(\d+)$',
        view_motif_reports.view_reports_by_TF_family),
]
