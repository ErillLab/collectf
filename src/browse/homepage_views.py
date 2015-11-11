from django.shortcuts import render

from . import motif_report


def home(request):
    """CollecTF homepage."""
    # Choose a random motif to be displayed on the main page.
    return render(request, 'homepage_home.html',
                  {'motif_report': motif_report.random_motif_report()})

def about(request):
    """Returns the 'about' page."""
    return render(request, 'homepage_about.html')


def browse(request):
    """Returns the 'browse' page."""
    return render(request, 'homepage_browse.html')


def search(request):
    """Returns the 'search' page."""
    return render(request, 'homepage_search.html')


def compare(request):
    """Returns the 'compare' page."""
    return render(request, 'homepage_compare.html')


def contribute(request):
    """Returns the 'contribute' page."""
    return render(request, 'homepage_contribute.html')


def stats(request):
    """Returns the 'stats' page."""
    return render(request, 'homepage_stats.html')


def links(request):
    """Returns the 'links' page."""
    return render(request, 'homepage_links.html')


def feedback(request):
    """Returns the 'feedback' page."""
    return render(request, 'homepage_feedback.html')


def cite(request):
    """Returns the 'cite' page."""
    return render(request, 'homepage_cite.html')


def acknowledgments(request):
    """Returns the 'acknowledgments' page."""
    return render(request, 'homepage_acknowledgments.html')
