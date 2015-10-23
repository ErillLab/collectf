from django.core.urlresolvers import reverse
from django.http import HttpRequest
from django.test import TestCase
from django.test.client import Client

class BrowseTFTest(TestCase):
    def setUp(self):
        self.client = Client()
                                 
    def test_browse_TF_view_success(self):
        response = self.client.get('/browse/browse_TF/')
        self.assertEqual(response.status_code, 200)

    def test_get_results_TF_family_URL_resolve_success(self):
        url = reverse('browse.browse_TF.get_results_TF_family', args=[1])
        self.assertEqual(url, '/browse/get_results_TF_family/1/')

    def test_get_results_TF_URL_resolve_success(self):
        url = reverse('browse.browse_TF.get_results_TF', args=[1])
        self.assertEqual(url, '/browse/get_results_TF/1/')
        
