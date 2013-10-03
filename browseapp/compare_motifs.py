from browse_base import *
from django import forms
from django.contrib.formtools.wizard.views import CookieWizardView
import search

class MotifComparisonForm(forms.Form):
    pass

FORMS = [("motif_a", MotifComparisonForm),
         ("motif_b", MotifComparisonForm)]

TEMPLATES = {"motif_a": "motif_compare_search.html",
             "motif_b": "motif_compare_search.html"}

class MotifComparisonWizard(CookieWizardView):
    def get_template_names(self):
        return [TEMPLATES[self.steps.current]]

    def get_context_data(self, form, **kwargs):
        # update context to include search form in the first two steps of motif
        # comparison
        context = super(MotifComparisonWizard, self).get_context_data(form=form, **kwargs)
        if self.steps.current in ["motif_a", "motif_b"]:
            template = search.search_get_template()
            context.update(template)
        return context

    def done(self, form_list, **kwargs):
        return render_to_response("success.html")


    
