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

    def test_get_results_taxonomy_URL_resolve_success(self):
        url = reverse('browse.browse_tax.get_results_taxonomy', args=[123])
        self.assertEqual(url, '/browse/get_results_taxonomy/123/')

