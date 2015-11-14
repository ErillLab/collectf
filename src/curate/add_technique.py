from django.shortcuts import redirect
from django.shortcuts import render
from django.contrib.auth.decorators import login_required

from .forms import add_technique as add_technique_form


@login_required
def add_technique(request):
    if request.method == 'POST':
        form = add_technique_form.ExperimentalTechniqueForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('homepage_home')

    else:
        form = add_technique_form.ExperimentalTechniqueForm()
    return render(request, 'add_technique.html', {'form': form})
