from django.core.urlresolvers import reverse
from django.http import HttpRequest
from django.test import TestCase
from django.test.client import Client

class MotifComparisonTest(TestCase):
    def setUp(self):
        self.client = Client()
                                 
    def test_compare_motif_page_step_1_view_success(self):
        response = self.client.get('/browse/compare_motifs/')
        self.assertEqual(response.status_code, 200)

    def test_compare_motif_page_step_2_view_success(self):
        response = self.client.get('/browse/compare_motifs/2/')
        self.assertEqual(response.status_code, 200)

    


        
