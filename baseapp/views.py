import bioutils
from django.http import HttpResponse

def get_weblogo(request):
    sites = request.POST['sites'].strip().split(',')
    sites = map(str, sites)
    assert all (len(sites[0]) == len(site) for site in sites)
    return HttpResponse(bioutils.weblogo_uri(sites), 'text/plain')
