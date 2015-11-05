from django.conf.urls import url

from .views import test_views

urlpatterns = [
    url(r'^test_view$', test_views.test_view),
]
