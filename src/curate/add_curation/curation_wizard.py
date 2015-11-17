"""Curation form wizard.

Here's the basic workflow for how a user would use a wizard:

1) The user visits the first page of the wizard, fills in the form and submits
it.

2) The server validates the data. If it's invalid, the form is displayed again,
with error messages. If it's valid, the server saves the current state of the
wizard in the backend and redirects to the next step.

3) Step 1 and 2 repeat, for every subsequent form in the wizard.

4) Once the user has submitted all the forms and all the data has been
validated, the wizard processes the data - saving it to the database, sending
an email, or whatever the application needs to do.

For more information, see django-formtools documentation.
"""

from django.contrib.auth.decorators import login_required
from formtools.wizard.views import SessionWizardView

from curate.forms.add_curation.publication import PublicationForm
from curate.forms.add_curation.genome import GenomeForm
from curate.forms.add_curation.techniques import TechniquesForm

from . import get_form
from . import process_step


class CurationWizard(SessionWizardView):
    """Form wizard to handle curation forms."""

    def get_template_names(self):
        """Overrides the default template names for each form."""
        return ['curation_step_%s.html' % self.steps.current]

    def get_form(self, step=None, data=None, files=None):
        """Constructs the form for a given step."""
        form = super(CurationWizard, self).get_form(step, data, files)
        if step is None:
            step = self.steps.current

        handlers = {
            '0': get_form.publication_get_form,
            '1': get_form.genome_get_form,
            '2': get_form.techniques_get_form,
        }
        return handlers[step](self, form)

    def process_step(self, form):
        """Post-processes data after each step."""
        handlers = {
            '0': process_step.publication_process_step,
            '1': process_step.genome_process_step,
            '2': process_step.techniques_process_step,
        }
        handlers[self.steps.current](self, form)
        return self.get_form_step_data(form)


@login_required
def curation(request):
    request.session['is_high_throughput_curation'] = False

    view = CurationWizard.as_view([
        PublicationForm,
        GenomeForm,
        TechniquesForm,
    ])

    return view(request)
