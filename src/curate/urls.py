from django.conf.urls import url

from . import add_TF
from . import add_technique

urlpatterns = [
    url(r'^add_TF/$', add_TF.add_TF, name='add_TF'),
    url(r'^add_TF_family/$', add_TF.add_TF_family, name='add_TF_family'),
    url(r'^add_technique/$', add_technique.add_technique,
        name='add_technique'),
]
