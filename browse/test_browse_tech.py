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

    def test_browse_by_function_invalid_function(self):
        assert False

    def test_browse_by_function_no_techniques(self):
        assert False

    def test_browse_by_function_no_curation_site_instances(self):
        assert False

    def test_browse_by_function_success(self):
        assert False

    def test_browse_by_category_invalid_function(self):
        assert False

    def test_browse_by_category_no_category(self):
        assert False

    def test_browse_by_category_no_techniques(self):
        assert False

    def test_browse_by_category_no_curation_site_instances(self):
        assert False

    def test_browse_by_category_success(self):
        assert False

    def test_browse_by_technique_no_technique(self):
        assert False

    def test_browse_by_technique_no_curation_site_instances(self):
        assert False

    def test_browse_by_tecnhique_success(self):
        assert False
        

    

