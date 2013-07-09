from browse_base import *

def browse_curation(request, cid):
    """Handler function for curation view"""
    curation = models.Curation.objects.get(curation_id=cid)
    return render_to_response("browse_curation.html",
                              {'curation': curation},
                              context_instance = RequestContext(request))


