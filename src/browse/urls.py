from django.conf.urls import url

from . import test_view
from . import homepage_views

urlpatterns = [
    url(r'^test_view$', test_view.test_view),
    url(r'^about$', homepage_views.about, name='homepage_about'),
    url(r'^browse$', homepage_views.browse, name='homepage_browse'),
    url(r'^search$', homepage_views.search, name='homepage_search'),
    url(r'^compare$', homepage_views.compare, name='homepage_compare'),
    url(r'^contribute$', homepage_views.contribute, name='homepage_contribute'),
    url(r'^stats$', homepage_views.stats, name='homepage_stats'),
    url(r'^links$', homepage_views.links, name='homepage_links'),
    url(r'^feedback$', homepage_views.feedback, name='homepage_feedback'),
    url(r'^cite$', homepage_views.cite, name='homepage_cite'),
    url(r'^acknowledgements$', homepage_views.acknowledgements,
        name='homepage_acknowledgements'),

    
]
