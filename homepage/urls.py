from django.conf.urls import patterns, url

urlpatterns = patterns('homepage.views',
    url(r'^greet/$', 'greet'),
    url(r'^about/$', 'about'),
    url(r'^browse/$', 'browse'),
    url(r'^search/$', 'search'),
    url(r'^compare/$', 'compare'),
    url(r'^contribute/$', 'contribute'),
    url(r'^feedback/$', 'feedback'),
    url(r'^stats/$', 'stats'),
    url(r'^cite/$', 'cite'),
    url(r'^links/$', 'links'),
    url(r'^acknowledgement/$', 'acknowledgements'),
    url(r'^feedback_send_email/$', 'feedback_send_email'),
    #url(r'^register/$', 'register_request'),
)

