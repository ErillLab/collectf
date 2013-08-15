import bioutils
from django.http import HttpResponse

def get_weblogo(request):
    sites = request.POST['sites'].strip().split(',')
    sites = map(str, sites)
    return HttpResponse(bioutils.weblogo_uri(sites), 'text/plain')
