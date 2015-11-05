from django.conf.urls import url

from . import test_view
from . import homepage_view

urlpatterns = [
    url(r'^test_view$', test_view.test_view),
    url(r'^home$', homepage_view.homepage),
]
