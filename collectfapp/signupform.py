from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

class CuratorRegistrationForm(UserCreationForm):
    first_name = forms.CharField(max_length=30)
    last_name = forms.CharField(max_length=30)
    email = forms.EmailField(max_length=75)

    class Meta:
        model = User
        fields = ("username", "first_name", "last_name", "email")

    def save(self, commit=True):
        cd = self.cleaned_data
        user = super(CuratorRegistrationForm, self).save(commit=False)
        user.first_name = cd["first_name"]
        user.last_name = cd["last_name"]
        user.email = cd["email"]
        if commit:
            user.save()
        return user
