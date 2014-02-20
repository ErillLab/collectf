"""
Handler class for curation process.

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

from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.contrib import messages
from django.contrib.formtools.wizard.views import SessionWizardView
from base import models # import models from base app
from forms.add_curation_form import * # import all curation forms
#from django.forms.formsets import formset_factory
import base
import session_utils
import add_curation_text
import add_curation_get_context_data
import add_curation_get_form
import add_curation_process_step
import add_curation_done

class CurationWizard(SessionWizardView):
    """Form wizard to handle curation forms. For methods and the session wizard
    logic, see Django docs:
    https://docs.djangoproject.com/en/1.6/ref/contrib/formtools/form-wizard/"""

    def get_template_names(self):
        """Override the default template name for each form"""
        return ["add_curation_step_%s.html" % self.steps.current]

    def get_context_data(self, form, **kwargs):
        """Return the template context for each step."""
        context = super(CurationWizard, self).get_context_data(form=form, **kwargs)
        context["form_title"] = add_curation_text.titles[self.steps.current]
        context["form_description"] = add_curation_text.descriptions[self.steps.current]
        # Call the associated function to update the context before rendering the step.
        funs = {'0': add_curation_get_context_data.publication_context_data,
                '1': add_curation_get_context_data.genome_context_data,
                '2': add_curation_get_context_data.techniques_context_data,
                '3': add_curation_get_context_data.site_entry_context_data,
                '4': add_curation_get_context_data.site_exact_match_context_data,
                '5': add_curation_get_context_data.site_soft_match_context_data,
                '6': add_curation_get_context_data.site_annotation_context_data,
                '7': add_curation_get_context_data.gene_regulation_context_data,
                '8': add_curation_get_context_data.curation_review_context_data,
        }
        context.update(funs[self.steps.current](self))
        return context

    def get_form(self, step=None, data=None, files=None):
        """Construct the form for a given step. The method is overriden to
        add/update arguments to the form instance."""
        form = super(CurationWizard, self).get_form(step, data, files)
        if step == None:
            step = self.steps.current
        
        handlers = {'0': add_curation_get_form.publication_get_form, # functions
                    '1': add_curation_get_form.genome_get_form,
                    '2': add_curation_get_form.techniques_get_form,
                    '3': add_curation_get_form.site_entry_get_form,
                    '4': add_curation_get_form.site_exact_match_get_form,
                    '5': add_curation_get_form.site_soft_match_get_form,
                    '6': add_curation_get_form.site_annotation_get_form,
                    '7': add_curation_get_form.gene_regulation_get_form,
                    '8': add_curation_get_form.curation_review_get_form,
        }
        return handlers[step](self, form)

    def process_step(self, form):
        """After each step, post-process data."""
        handlers = {'0': add_curation_process_step.publication_process, # functions
                    '1': add_curation_process_step.genome_process,
                    '2': add_curation_process_step.techniques_process,
                    '3': add_curation_process_step.site_entry_process,
                    '4': add_curation_process_step.site_exact_match_process,
                    '5': add_curation_process_step.site_soft_match_process,
                    '6': add_curation_process_step.site_annotation_match_process,
                    '7': add_curation_process_step.gene_regulation_process,
                    '8': add_curation_process_step.review_curation_process,
        }
        handlers[self.steps.current](self, form)
        return self.get_form_step_data(form)

    def render(self, form=None, **kwargs):
        """
        This method gets called after the GET or POST request has been handled.

        The only reason to override this method is to check "no data" button in
        publication selection form step. If the curator selects a publication
        and checks the button "This paper contains no data.", the proper action
        is to terminate curation process as ther is no data to be curated. The
        problem is to process that information and redirect to the home
        page. This is achieved by overloading render method. Before this
        function is called, publication object is modified as having no
        data. Afterwards, this function is called and it redirects to the
        homepage with the message about the action that was performed.
        """
        
        form = form or self.get_form()
        context = self.get_context_data(form=form, **kwargs)

        if (session_utils.has(self.request.session, "paper_contains_no_data") and
            session_utils.get(self.request.session, "paper_contains_no_data")):
            msg = "The publication was marked as having no data."
            messages.info(self.request, msg)
            # clear session data
            session_utils.clear(self.request.session)
            return HttpResponseRedirect(reverse(base.views.home))
        
        return self.render_to_response(context)

    def done(self, form_list, **kwargs):
        """Last step in the curation process. This method is called after all
        forms. Insert all data into the database."""
        return add_curation_done.master_done(self, form_list, **kwargs)

@login_required
def curation(request):
    view = CurationWizard.as_view([PublicationForm,
                                   GenomeForm,
                                   TechniquesForm,
                                   SiteEntryForm,
                                   SiteExactMatchForm,
                                   SiteSoftMatchForm,
                                   SiteAnnotationForm,
                                   GeneRegulationForm,
                                   CurationReviewForm,],
                                  condition_dict = {'5': inexact_match_form_condition})
    return view(request)

def inexact_match_form_condition(wizard):
    """Check if inexact match form is necessary. If not (i.e. all sites have
    been matched exactly, hide this step.)"""
    sites = session_utils.get(wizard.request.session, 'sites')
    return any(not site.is_matched() for site in sites)

