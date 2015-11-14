from django.forms import ModelForm

from core import models


class TFForm(ModelForm):
    class Meta:
        model = models.TF
        fields = '__all__'


class TFFamilyForm(ModelForm):
    class Meta:
        model = models.TFFamily
        fields = '__all__'
