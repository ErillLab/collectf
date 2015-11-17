from django import forms
from . import help_text

class PublicationForm(forms.Form):
    """Publication selection form.

    The curators are asked to select one of the papers assigned to them.
    """
    help_text_dict = help_text.publication_form
    
    publication = forms.ChoiceField(
        label="Publications",
        widget=forms.RadioSelect())

    no_data = forms.BooleanField(
        label="This paper contains no data.",
        required=False,
        help_text=help_text_dict['no_data'])
