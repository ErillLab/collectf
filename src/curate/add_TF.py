from django.shortcuts import redirect
from django.shortcuts import render
from django.contrib.auth.decorators import login_required

from .forms import add_TF as add_TF_form


@login_required
def add_TF(request):
    if request.method == 'POST':
        form = add_TF_form.TFForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('homepage_home')

    else:
        form  = add_TF_form.TFForm()
    return render(request, 'add_TF.html', {'form': form})


@login_required
def add_TF_family(request):
    if request.method == 'POST':
        form = add_TF_form.TFFamilyForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('homepage_home')

    else:
        form  = add_TF_form.TFFamilyForm()
    return render(request, 'add_TF_family.html', {'form': form})
