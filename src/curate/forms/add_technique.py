from django.forms import ModelForm

from core import models


class ExperimentalTechniqueForm(ModelForm):
    class Meta:
        model = models.ExperimentalTechnique
        fields = '__all__'
