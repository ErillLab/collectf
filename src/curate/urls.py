from django.conf.urls import url

from . import add_TF
from . import add_publication
from . import add_technique
from . import json_views
from . import add_curation

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

    # Curation
    url(r'^curation/$', add_curation.curation, name='curation'),
    url(r'^high_throughput_curation/$', add_curation.high_throughput_curation,
        name='high_throughput_curation'),

    # Helper JSON views
    url(r'^json_get_genomes$', json_views.get_genomes),
    url(r'^json_get_TF_instances$', json_views.get_TF_instances),
    url(r'^json_uniprot_to_refseq/(\w+)$', json_views.uniprot_to_refseq),
]
