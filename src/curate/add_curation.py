"""
Handler class for curation process.

From Django docs:

Here's the basic workflow for how a user would use a wizard:

1. The user visits the first page of the wizard, fills in the form and submits
it.

2. The server validates the data. If it's invalid, the form is displayed again,
with error messages. If it's valid, the server saves the current state of the
wizard in the backend and redirects to the next step.

3. Step 1 and 2 repeat, for every subsequent form in the wizard.

4. Once the user has submitted all the forms and all the data has been
validated, the wizard processes the data - saving it to the database, sending an
email, or whatever the application needs to do.
"""

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from formtools.wizard.views import SessionWizardView
from django.core.urlresolvers import reverse
from django.shortcuts import redirect

from curate import add_curation_done
from curate import add_curation_get_context_data
from curate import add_curation_get_form
from curate import add_curation_process_step
from curate.forms.add_curation_form import CurationReviewForm
from curate.forms.add_curation_form import GeneRegulationForm
from curate.forms.add_curation_form import GenomeForm
from curate.forms.add_curation_form import PublicationForm
from curate.forms.add_curation_form import SiteAnnotationForm
from curate.forms.add_curation_form import SiteEntryForm
from curate.forms.add_curation_form import SiteExactMatchForm
from curate.forms.add_curation_form import SiteSoftMatchForm
from curate.forms.add_curation_form import TechniquesForm

class CurationWizard(SessionWizardView):
    """Form wizard to handle curation forms.

    For methods and the session wizard logic, see Django docs:
    https://docs.djangoproject.com/en/1.6/ref/contrib/formtools/form-wizard/
    """

    def get_template_names(self):
        """Override the default template name for each form"""
        return ["curation_step_%s.html" % self.steps.current]

    def get_context_data(self, form, **kwargs):
        """Return the template context for each step."""
        context = super(CurationWizard, self).get_context_data(form=form,
                                                               **kwargs)
        # Call the associated function to update the context before rendering
        # the step.
        funs = {
            '0': add_curation_get_context_data.publication_context_data,
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
        """Constructs the form for a given step.

        The method is overridden to add/update arguments to the form
        instance.
        """
        form = super(CurationWizard, self).get_form(step, data, files)
        if step == None:
            step = self.steps.current

        handlers = {'0': add_curation_get_form.publication_get_form,
                    '1': add_curation_get_form.genome_get_form,
                    '2': add_curation_get_form.techniques_get_form,
                    '3': add_curation_get_form.site_entry_get_form,
                    '4': add_curation_get_form.site_exact_match_get_form,
                    '5': add_curation_get_form.site_soft_match_get_form,
                    '6': add_curation_get_form.site_annotation_get_form,
                    '7': add_curation_get_form.gene_regulation_get_form,
                    '8': add_curation_get_form.curation_review_get_form}
        return handlers[step](self, form)

    def process_step(self, form):
        """Post-processes data after each step."""
        handlers = {'0': add_curation_process_step.publication_process,
                    '1': add_curation_process_step.genome_process,
                    '2': add_curation_process_step.techniques_process,
                    '3': add_curation_process_step.site_entry_process,
                    '4': add_curation_process_step.site_exact_match_process,
                    '5': add_curation_process_step.site_soft_match_process,
                    '6': add_curation_process_step.site_annotation_process,
                    '7': add_curation_process_step.gene_regulation_process,
                    '8': add_curation_process_step.review_curation_process}
        handlers[self.steps.current](self, form)
        return self.get_form_step_data(form)

    def render(self, form=None, **kwargs):
        """Gets called after the GET or POST request has been handled.

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

        session = self.request.session
        if (session.get('paper_contains_no_data') and
            session.get('paper_contains_no_data')):
            msg = "The publication was marked as having no data."
            messages.info(self.request, msg)
            session.clear() # clear session data
            return redirect('homepage_home')

        return self.render_to_response(context)

    def done(self, form_list, **kwargs):
        """Gets called after all forms. Inserts all data into the database."""
        return add_curation_done.master_done(self, form_list, **kwargs)

@login_required
def curation(request):
    """Entry point for the curation."""

    # If user selects the old curation and then go back, the session will have
    # the old_curation key in table, and it will cause trouble.
    if request.session.get('old_curation'):
        request.session['old_curation'] = None

    # This is not high-throughput submission.
    request.session['high_throughput_curation'] = False

    view = CurationWizard.as_view(
        [PublicationForm,
         GenomeForm,
         TechniquesForm,
         SiteEntryForm,
         SiteExactMatchForm,
         SiteSoftMatchForm,
         SiteAnnotationForm,
         GeneRegulationForm,
         CurationReviewForm,],
        condition_dict={'5': inexact_match_form_condition})
    return view(request)

@login_required
def high_throughput_curation(request):
    """Entry point for high-throughput curation.

    Curators can check ChIP and other high-throughput methodologies in the
    regular submission mode, but if they are submitting data that is primarily
    based on high-throughput binding assays (e.g. ChIP-seq genomic-SELEX, etc.)
    they are then encouraged to use the high-throughput submission portal. First
    few steps are identical to the ones in the regular submission portal. In the
    site-entry step, two types of data are asked: sites and peaks. As in regular
    submission portal, sites can be either motif-associated,
    non-motif-associated or variable-motif-associated (e.g. variable spacing,
    inverting, anything that is not gapless alignment), which can be entered in
    sequence-based or coordinate-based modes. Below the site box, curators are
    able to enter peak data (most likely in coordinate mode).
    """

    # This IS high-throughput submission
    request.session['high_throughput_curation'] = True

    view = CurationWizard.as_view(
        [PublicationForm,
         GenomeForm,
         TechniquesForm,
         SiteEntryForm,
         SiteExactMatchForm,
         SiteSoftMatchForm,
         SiteAnnotationForm,
         GeneRegulationForm,
         CurationReviewForm,],
        condition_dict={'5': inexact_match_form_condition})
    return view(request)

def inexact_match_form_condition(wizard):
    """Checks if inexact match form is necessary.

    If not, it means that all sites have been matched exactly; hide this step.
    """
    sites = wizard.request.session.get('sites')
    all_exact_matched = (sites and
                         all(site.is_matched() and site.get_match().is_exact()
                             for site in sites))
    return not all_exact_matched
