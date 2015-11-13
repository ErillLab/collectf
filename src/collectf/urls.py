"""collectf URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/dev/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Import the include() function: from django.conf.urls import url, include
    3. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
from django.conf.urls import url, include
from django.contrib import admin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import RedirectView

from browse import urls as browse_urls
from browse import view_site

urlpatterns = [
    url(r'^$', RedirectView.as_view(url=reverse_lazy('homepage_home'))),
    url(r'^admin/', admin.site.urls),
    url(r'^browse/', include(browse_urls)),

    # dbxref links from NCBI don't have browse/ prefix, being served from here.
    url(r'^expsite_(?P<dbxref_id>\w+)$', view_site.view_site),
    url(r'^EXPSITE_(?P<dbxref_id>\w+)$', view_site.view_site),

    # user account management
    url(r'^accounts/', include('registration.backends.default.urls')),
]
