from django.conf.urls import url

from . import add_TF
from . import add_technique
from . import add_publication

urlpatterns = [
    # Add TF
    url(r'^add_TF/$', add_TF.add_TF, name='add_TF'),
    
    # Add TF family
    url(r'^add_TF_family/$', add_TF.add_TF_family, name='add_TF_family'),

    # Add experimental technique
    url(r'^add_technique/$', add_technique.add_technique,
        name='add_technique'),

    # Add PubMed publication
    url(r'^add_pubmed_publication/$', add_publication.add_pubmed_publication,
        name='add_pubmed_publication'),
    # Add non-PubMed publication
    url(r'^add_non_pubmed_publication/$',
        add_publication.add_non_pubmed_publication,
        name='add_non_pubmed_publication'),
]
