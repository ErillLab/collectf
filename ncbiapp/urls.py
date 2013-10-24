from django.conf.urls import patterns, url

urlpatterns = patterns('ncbiapp.views',
    url(r'^export_tbl_view', 'export_tbl_view'),
)


