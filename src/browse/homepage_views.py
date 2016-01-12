from django.contrib import messages
from django.core.mail import send_mail
from django.shortcuts import render
from django.shortcuts import redirect

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


def feedback_send_email(request):
    """Sends an email to CollecTF team containing a feedback form."""
    try:
        send_mail('CollecTF feedback (%s)' % request.POST['type'],
                  request.POST['comment'],
                  request.POST['email'],
                  ['sefa1@umbc.edu', 'collectfdb@umbc.edu'],
                  fail_silently=False)
        messages.add_message(request, messages.INFO,
                             "Thanks for the feedback!")
    except:
        messages.add_message(
            request, messages.ERROR,
            """Something went wrong. Please try again or send your feedback
            directly to collectfdb@umbc.edu""")

    return redirect('homepage_home')


def cite(request):
    """Returns the 'cite' page."""
    return render(request, 'homepage_cite.html')


def acknowledgments(request):
    """Returns the 'acknowledgments' page."""
    return render(request, 'homepage_acknowledgments.html')
