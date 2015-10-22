from django.contrib.auth.models import AnonymousUser
from django.core.urlresolvers import resolve
from django.http import HttpRequest
from django.test import TestCase

import views

class HomePageTest(TestCase):
    def setUp(self):
        self.request = HttpRequest()
        self.request.user = AnonymousUser()
        
    def test_root_url_resolves_to_home_page_view(self):
        found = resolve('/')
        self.assertEqual(found.func, views.home)
