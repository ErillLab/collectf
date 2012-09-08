"""Handler class for curation process.

From Django docs:

Here's the basic workflow for how a user would use a wizard:

1) The user visits the first page of the wizard, fills in the form and submits it.
2) The server validates the data. If it's invalid, the form is displayed again,
with error messages. If it's valid, the server saves the current state of the
wizard in the backend and redirects to the next step.
3) Step 1  and 2 repeat, for every subsequent form in the wizard.
4) Once the user has submitted all the forms and all the data has been
validated, the wizard processes the data - saving it to the database, sending an
email, or whatever the application needs to do.
"""

from django.contrib.formtools.wizard.views import SessionWizardView
from curationform import *
import models
import sutils

# curation get form functions
# get_form constructs the form for a given step
def publication_form(wiz, form):
    """Publication selection form"""
    user = wiz.request.user
    curator = models.Curator.objects.get(user=user)
    # select papers which are not complete yet
    assigned_pubs = curator.assigned_papers.filter(curation_complete=False)
    # put them in form choices, populate form field
    choices = [(p.publication_id, p.citation) for p in assigned_pubs]
    form.fields["pub"].choices = choices

    return form
    
def genome_form(wiz, form):
    """Genome, TF, TF_family, tf instance, .. selection form"""
    # populate TF field
    choices = [(None, 'None'),]
    for TF in models.TF.objects.all():
        choices.append((TF.TF_id, TF.name + " (family: %s)" % TF.family))
    form.fields['TF'].choices = choices
    return form

# curation process step functions
# Post process the form data before the data gets stored
def publication_process(wiz, form):
    pubid = form.cleaned_data['pub']
    sutils.sput(wiz.request.session, 'publication', pubid)

def genome_process(wiz, form):
    genome_accession = form.cleaned_data['genome_accession']
    # in form validation genome is searched in db, and if not found,
    # it is inserted into db. So, at this point, it is guaranteed that
    # genome with id <genome_accession> should be in db.
    genome = models.Genome.objects.get(genome_accession=genome_accession)

class CurationWizard(SessionWizardView):
    """Form wizard to handle curation forms. For methods, see Django docs."""

    def get_template_names(self):
        """Override default template name for each form, default is
        wizard_form.html"""
        return ['wizard_form.html']

    def get_form(self, step=None, data=None, files=None):
        """Construct the form for a given step. The method is overriden to
        add/update arguments to the form instance."""
        form = super(CurationWizard, self).get_form(step, data, files)
        if step == None:
            step = self.steps.current

        handlers = {'0': publication_form, # functions
                    '1': genome_form,
                   }
        
        return handlers[step](wiz=self, form=form)

    def process_step(self, form):
        """Process data after each step before it gets stored"""
        handlers = {'0': publication_process, # functions
                    '1': genome_process,
                   }
        
        handlers[self.steps.current](wiz=self, form=form)
        return self.get_form_step_data(form)
        


    
# curation handler
# for form definitions, go curationform.py
curation = CurationWizard.as_view([PublicationForm,
                                   GenomeForm,])
