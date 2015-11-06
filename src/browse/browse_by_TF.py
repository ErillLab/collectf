"""The view functions for browsing by TF and TF family."""

from django.shortcuts import render

from core import models

def browse_by_TF(request):
    """Returns the TF treeview for browsing by TF and family."""
    context = {'TF_families': models.TFFamily.objects.all().order_by('name')}
    print context
    return render(request, 'browse_by_TF.html', context)

