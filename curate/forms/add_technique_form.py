"""Forms for adding experimental techniques"""

from django import forms
import base.models as models

class AddTechniqueForm(forms.Form):
    name = forms.CharField(label="Name")
    description = forms.CharField(label="Description", widget=forms.Textarea)
    type = forms.ChoiceField(choices = models.ExperimentalTechnique.FUNCTION_CATEGORIES)
    categories = forms.ModelMultipleChoiceField(queryset=models.ExperimentalTechniqueCategory.objects.order_by('name'))
