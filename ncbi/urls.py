from django.conf.urls import patterns, url

urlpatterns = patterns('ncbi.views',
    url(r'^export_tbl_view', 'export_tbl_view'),
)

