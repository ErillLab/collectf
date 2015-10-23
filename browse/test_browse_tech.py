from django.core.urlresolvers import reverse
from django.http import HttpRequest
from django.test import TestCase
from django.test.client import Client

class BrowseTechniquesTest(TestCase):

    def setUp(self):
        self.client = Client()
        
    def test_browse_techniques_view_success(self):
        response = self.client.get('/browse/browse_technique/')
        self.assertEqual(response.status_code, 200)

    def test_get_results_technique_URL_resolve_success(self):
        url = reverse('browse.browse_tech.get_results_technique', args=[123])
        self.assertEqual(url, '/browse/get_results_technique/123/')

    def test_get_results_technique_all_URL_resolve_success(self):
        url = reverse('browse.browse_tech.get_results_all', args=['binding'])
        self.assertEqual(url, '/browse/get_results_technique_all/binding/')

    def test_get_results_technique_category_URL_resolve_success(self):
        url = reverse('browse.browse_tech.get_results_category',
                      args=['binding', 1])
        self.assertEqual(
            url, '/browse/get_results_technique_category/binding/1/')

    

