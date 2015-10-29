from django.core.urlresolvers import reverse
from django.http import HttpRequest
from django.test import TestCase
from django.test.client import Client

class BrowseTaxTest(TestCase):

    def setUp(self):
        self.client = Client()
        
    def test_browse_taxonomy_view_success(self):
        response = self.client.get('/browse/browse_taxonomy/')
        self.assertEqual(response.status_code, 200)

    def test_browse_taxonomy_no_taxonomy(self):
        assert False

    def test_browse_taxonomy_empty_curation_site_instances(self):
        assert False

    def test_browse_taxonomy_success(self):
        assert False
