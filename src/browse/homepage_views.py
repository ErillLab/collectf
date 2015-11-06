from django.shortcuts import render


def home(request):
    """CollecTF homepage."""
    return render(request, 'homepage_home.html')


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


def acknowledgements(request):
    """Returns the 'acknowledgements' page."""
    return render(request, 'homepage_acknowledgements.html')
