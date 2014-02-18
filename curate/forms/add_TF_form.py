"""Forms for adding a TF or TF family"""

from django import forms
import base.models as models

class AddTFForm(forms.Form):
    name = forms.CharField(label="Name")
    family = forms.ModelChoiceField(queryset=models.TFFamily.objects.order_by('name'))
    description = forms.CharField(widget=forms.Textarea)

class AddTFFamilyForm(forms.Form):
    name = forms.CharField(label="Name")
    description = forms.CharField(widget=forms.Textarea)
