from django.core.urlresolvers import reverse
from django.http import HttpRequest
from django.test import TestCase
from django.test.client import Client

from . import browse_TF
from base import models

class BrowseTFTest(TestCase):
    def setUp(self):
        self.client = Client()
        self.TF_family = models.TFFamily.objects.create(name='TF_family')
        self.TF = models.TF.objects.create(name='TF', family=TF_family)
                                 
    def test_browse_TF_view_success(self):
        response = self.client.get('/browse/browse_TF/')
        self.assertEqual(response.status_code, 200)

    def test_browse_TF_family_no_TF_family(self):
        assert False

    def test_browse_TF_family_no_TF(self):
        assert False
    
    def test_browse_TF_family_no_curation_site_instances(self):
        assert False

    def test_browse_TF_family_success(self):
        assert False

    def test_browse_TF_no_TF(self):
        assert False

    def test_browse_TF_no_curation_site_instances(self):
        assert False

    def test_browse_TF_success(self):
        assert False
