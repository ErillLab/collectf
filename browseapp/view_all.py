# Handler for viewing all curations in a page
# Handler for viewing all publications in a page

from browse_base import *


@login_required
def view_all_curations(request):
    """Handler function to see all curations at once.
    This function renders the page with the list of
    all curations in the database"""
    all_curations = models.Curation.objects.all().order_by('curation_id')
    return render_to_response("view_all_curation.html",
                              {"curations": all_curations},
                              context_instance=RequestContext(request))

@login_required
def view_all_publications(request):
    """Handler function to see all publications in the database.
    This is for internal use, to see all publications in the database."""
    all_pubs = models.Publication.objects.all().order_by('-pmid')
    return render_to_response("view_all_publication.html",
                              {"publications": all_pubs},
                              context_instance=RequestContext(request))
