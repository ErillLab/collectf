from django.conf.urls import patterns, url


urlpatterns = patterns('mastercuratorapp.views',
    url(r'^home/$', 'home'),
)
