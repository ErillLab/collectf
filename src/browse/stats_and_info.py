"""View functions about CollecTF stats and info."""

from django.shortcuts import render

def curator_roster(request):
    """Returns the list of curators."""
    return render(request, 'curator_roster.html')
