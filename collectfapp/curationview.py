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

from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.contrib.formtools.wizard.views import SessionWizardView
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from curationform import *
import models
import sutils
import sitesearch
import views

# curation get form functions
# get_form constructs the form for a given step
def publication_get_form(wiz, form):
    """Publication selection form"""
    user = wiz.request.user
    curator = models.Curator.objects.get(user=user)
    # select papers assigned to user
    assigned_pubs = models.Publication.objects.filter(assigned_to=curator)
    # select papers which are not complete yet
    not_completed_pubs = assigned_pubs.filter(curation_complete=False)
    # put them in form choices, populate form field
    choices = [(p.publication_id, "[%s] %s" % (p.publication_id, p.citation))
                for p in not_completed_pubs]
    form.fields["pub"].choices = choices
    return form
    
def genome_get_form(wiz, form):
    """Genome, TF, TF_family, tf instance, .. selection form"""
    # populate TF field
    return form

def techniques_get_form(wiz, form):
    # populate the following fields
    # - contains expression data
    # - contains promoter data
    pid = sutils.sget(wiz.request.session, 'publication')
    pub = models.Publication.objects.get(publication_id=pid)
    form.fields["contains_promoter_data"].initial = pub.contains_promoter_data
    form.fields["contains_expression_data"].initial = pub.contains_expression_data
    return form

def site_report_get_form(wiz, form):
    return form

# helper function for site_exact_match_form and site_soft_match_forms
def populate_match_choices(sid, matches):
    # for a given site and its all matches, populate Django field choice
    # exact_match: True if populating exact match results, if false, soft search
    choices = []
    for mid, match in matches.items():
        choices.append((mid, sitesearch.print_site_match(match)))
    last_choice_msg = "No valid match"
    choices.append((None, last_choice_msg))
    return choices
    
def site_exact_match_get_form(wiz, form):
    """Show list of sites and their exact matches."""
    # Each reported site and its exact match results are represented as field.
    # Dynamically add fields to site search results step.
    sites = sutils.sget(wiz.request.session, 'sites')
    site_match_choices = sutils.sget(wiz.request.session, 'site_match_choices')
    for sid, matches in site_match_choices.items(): # for all matches belong to a site
        choices = populate_match_choices(sid, matches)
        # make the form field
        form.fields[sid] = forms.ChoiceField(label=sites[sid], choices=choices,
                                             widget=forms.RadioSelect())
    return form

def site_soft_match_get_form(wiz, form):
    """For sites that are not matched with equivalents in the previous form,
    show soft search results."""
    sites = sutils.sget(wiz.request.session, 'sites')
    # soft site matches: {sid: {mid: SiteMatch}}
    soft_site_match_choices = sutils.sget(wiz.request.session, 'soft_site_match_choices')
    # exact site matches: {sid: SiteMatch} -- they're already matched in prev form
    exact_site_matches = sutils.sget(wiz.request.session, 'exact_site_matches')
    for sid, matches in soft_site_match_choices.items():
        choices = populate_match_choices(sid, matches)
        # make the form field
        form.fields[sid] = forms.ChoiceField(label=sites[sid], choices=choices,
                                             widget=forms.RadioSelect())
    return form

def site_regulation_get_form(wiz, form):
    """For each site and its equivalent, show nearby genes for regulation input. For
    each nearby gene for each site, curator manually checks if site regulates gene or
    not (i.e. if there is any experimantal support for regulation)."""
    pubid = sutils.sget(wiz.request.session, 'publication')
    publication = models.Publication.objects.get(publication_id=pubid)  # get pub

    # get exact and soft site matches
    exact_site_matches = sutils.sget(wiz.request.session, 'exact_site_matches')
    soft_site_matches = sutils.sget(wiz.request.session, 'soft_site_matches')

    for sid,match in exact_site_matches.items() + soft_site_matches.items():
        choices = []  # list of genes to the site match
        for g in match.nearby_genes:
            choices.append((g.gene_id, '%s locus_tag: %s' % (g.name, g.locus_tag)))
        form.fields[sid] = forms.MultipleChoiceField(label=match.match.seq,
                                                     choices=choices, required=False,
                                                     widget=forms.CheckboxSelectMultiple)
        # disable checkbox if publication is marked as not having expression data
        if not publication.contains_expression_data:
            form.fields[sid].widget.attrs['disabled'] = 'disabled'
    return form

def curation_review_get_form(wiz, form):
    """Get curation review form"""
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
    # store genome in session data
    sutils.sput(wiz.request.session, 'genome', genome)

def techniques_process(wiz, form):
    # set publication fields
    # - contains_promoter_data
    # - contains_expression_data
    pubid = sutils.sget(wiz.request.session, 'publication')
    p = models.Publication.objects.get(publication_id=pubid)
    p.contains_promoter_data = form.cleaned_data["contains_promoter_data"]
    p.contains_expression_data = form.cleaned_data["contains_expression_data"]
    p.save()


def site_report_process(wiz, form):
    # read genome from session data
    genome = sutils.sget(wiz.request.session, 'genome')
    genes = models.Gene.objects.filter(genome=genome).order_by('start')
    # get reported sites -- list of (sid, site)
    sites = sitesearch.parse_site_input(form.cleaned_data['sites'])
    # find exact matches
    site_match_choices = sitesearch.match_all_exact(genome, genes, sites)
    # put sites and site matches into session data
    sutils.sput(wiz.request.session, 'site_match_choices', site_match_choices)
    sutils.sput(wiz.request.session, 'sites', sites)
    
def site_exact_match_process(wiz, form):
    """In the form, reported sites and their matches are displayed. For some
    sites, curator may choose to perform a 'soft' search (if surrogate genome is
    being used, sites may not be exactly found in the genome)."""
    genome = sutils.sget(wiz.request.session, 'genome')      # genome
    genes = models.Gene.objects.filter(genome=genome).order_by('start')
    sites = sutils.sget(wiz.request.session, 'sites')        # all sites
    # get all match choices for each reported site
    site_match_choices = sutils.sget(wiz.request.session, 'site_match_choices')
    exact_site_matches = {} # {sid: SiteMatch} to store which sites are exactly matched
    for sid, mid in form.cleaned_data.items():  # siteid and matchid
        if mid != 'None':
            exact_site_matches[sid] = site_match_choices[sid][int(mid)]
    # the rest is unmatched, perform soft search for them
    soft_sites = {}
    for sid, site in sites.items():
        if sid not in exact_site_matches:
            soft_sites[sid] = site

    exact_sites = [m.match.seq for m in exact_site_matches.values()]
    soft_site_match_choices = sitesearch.match_all_soft(genome, genes,
                                                        soft_sites, exact_sites)
    
    # - exact_site_matches {sid: SiteMatch}
    # - soft_site_match_choices {sid: {mid: SiteMatch}} and get user selections for those
    sutils.sput(wiz.request.session, 'exact_site_matches', exact_site_matches)
    sutils.sput(wiz.request.session, 'soft_site_match_choices', soft_site_match_choices)

def site_soft_match_process(wiz, form):
    """In this form, soft search results are processed, user matched all of sites to
    any equivalent. Some sites might be unmatched if there is not any good candidate
    in search results."""
    sites = sutils.sget(wiz.request.session, 'sites')
    soft_site_match_choices = sutils.sget(wiz.request.session, 'soft_site_match_choices')
    soft_site_matches = {} # {sid: SiteMatch} to store which sites are softly matched
    not_matched_sites = [] # list of ids of sites that are not matched
    for sid, mid in form.cleaned_data.items(): # siteid and matchid
        if mid != 'None':
            # reported site mapped to a soft search result site
            soft_site_matches[sid] = soft_site_match_choices[sid][int(mid)]
        else:
            # site is not mapped as exact nor soft search, insert to DB as unmatched
            not_matched_sites.append(sid)

    # - soft_site_matches {sid: SiteMatch}
    # - not_matched_sites [sid, ..]
    sutils.sput(wiz.request.session, 'soft_site_matches', soft_site_matches)
    sutils.sput(wiz.request.session, 'not_matched_sites', not_matched_sites)

    

def site_regulation_process(wiz, form):
    pass

def curation_review_process(wiz, form):
    pass

# CurationWizard done step functions

# When the last step of the form is submitted, gather all data from each form
# step and insert into database. Following _done_ functions are implemented for
# each step and they are called by CurationWizard.done when user submit curation
# form.
def publication_done(wiz, form, **kwargs):
    pid = form.cleaned_data['pub']
    publication = models.Publication.objects.get(publication_id=pid)
    return publication

def genome_done(wiz, form, **kwargs):
    """Get genome and TF data from the genome form (2nd step)"""
    cd = {}  # cleaned data to be returned
    # genome
    genome_accession = form.cleaned_data['genome_accession']
    cd['genome'] = models.Genome.objects.get(genome_accession=genome_accession)
    # TF instance
    TF_accession = form.cleaned_data['TF_accession']
    cd['TF_instance'] = models.TFInstance.objects.get(protein_accession=TF_accession)
    # get all other cleaned data
    cd['TF'] = form.cleaned_data['TF']
    cd['TF_function'] = form.cleaned_data['TF_function']
    cd['TF_type'] = form.cleaned_data['TF_type']
    cd['TF_species'] = form.cleaned_data['TF_species']
    cd['site_species'] = form.cleaned_data['site_species']
    return cd

def techniques_done(wiz, form, **kwargs):
    """Get techniques and experimental process data from form"""
    cd = {}  # cleaned data to be returned
    cd['techniques'] = form.cleaned_data['techniques']
    cd['experimental_process'] = form.cleaned_data['experimental_process']
    cd['forms_complex'] = form.cleaned_data['forms_complex']
    cd['complex_notes'] = form.cleaned_data['complex_notes']
    return cd

def site_report_done(wiz, form, **kwargs):
    pass

def site_exact_match_done(wiz, form, **kwargs):
    """Get data from exact site match form"""
    pass

def site_soft_match_done(wiz, form, **kwargs):
    """Get data from soft site match form"""
    pass

def site_regulation_done(wiz, form, **kwargs):
    """Get regulation data"""
    return form.cleaned_data

def curation_review_done(wiz, form, **kwargs):
    """Get review form data and return it"""
    cd = {}  # cleaned data to be returned
    cd['requires_revision'] = form.cleaned_data['revision_reasons']
    cd['confidence'] = form.cleaned_data['confidence']
    cd['NCBI_submission_ready'] = form.cleaned_data['NCBI_submission_ready']
    cd['paper_complete'] = form.cleaned_data['paper_complete']
    cd['notes'] = form.cleaned_data['notes']
    return cd

# other done helper functions

def regulation_done(wiz, regulations, match, curation_site_instance):
    """For site instance, add regulation information if available"""
    for gene in match.nearby_genes:
        ev = None
        if str(gene.gene_id) in regulations:
            ev = 'exp_verified'  # experimentally verified
        else:
            ev = 'inferred'  # inferred
        # add regulation object
        r = models.Regulation(curation_site_instance=curation_site_instance,
                              gene=gene, evidence_type=ev)
        r.save()

def site_match_done(wiz, curation, regulations):
    """Get exact and soft site matches and create those objects in the database."""
    sites = sutils.sget(wiz.request.session, 'sites')
    exact_site_matches = sutils.sget(wiz.request.session, 'exact_site_matches')
    soft_site_matches = sutils.sget(wiz.request.session, 'soft_site_matches')
    genome = sutils.sget(wiz.request.session, 'genome')  # get genome
    # combine exact and soft site matches
    all_site_matches = dict(exact_site_matches.items() + soft_site_matches.items())
    for sid, match in all_site_matches.items():
        # create SiteInstance object (or get if available)
        si, created  = models.SiteInstance.objects.get_or_create(
            genome=genome, start=match.match.start, end=match.match.end,
            strand=match.match.strand, seq=sites[sid])
                
        # create curation_siteinstance object (through relation)
        cs = models.Curation_SiteInstance(curation=curation, site_instance=si,
                                          annotated_seq=match.match.seq)
        cs.save()

        # for site instance, add regulation information
        regulation_done(wiz, regulations=regulations[sid], match=match,
                        curation_site_instance=cs)


        
def not_matched_sites_done(wiz, curation):
    """Create not matched site instances."""
    sites = sutils.sget(wiz.request.session, 'sites')
    not_matched_sites = sutils.sget(wiz.request.session, 'not_matched_sites')
    for sid in not_matched_sites:
        models.NotAnnotatedSiteInstance(sequence=sites[sid], curation=curation).save()
  

class CurationWizard(SessionWizardView):
    """Form wizard to handle curation forms. For methods, see Django docs."""

    def get_template_names(self):
        """Override default template name for each form, default is
        wizard_form.html"""
        return ['wizard_form.html']

    def get_context_data(self, form, **kwargs):
        """Return template context for a step."""
        context = super(CurationWizard, self).get_context_data(form=form, **kwargs)
        titles = {
            '0': "Publication selection",
            '1': "Genome and TF information",
            '2': "Experimental methods used in this paper",
            '3': "Reported sites",
            '4': "Exact site matches",
            '5': "Inexact site matches",
            '6': "Gene regulation (experimental support)",
            '7': "Curation information"
        }
        descriptions = {
            '0': "Please choose a publication to curate.",
            '1': """This step collects information on the transcription factor (TF), the specific
strains reported in the manuscript and the NCBI GenBank sequences that reported
sites and TF will be mapped onto.""",
            '2': """Select the experimental techniques and the describe the basic experimental
procedure used to verify binding/expression of the sites reported in this
curation.""",
            '3': """Enter the list of sites as reported in the paper. Supported formats are:\n(a)
Raw sequence (one site per line)\n(b) FASTA format""",
            '4': """For each reported site, all exact matches in the chosen genome are listed. If a
reported site does not have any exact matches, or the matched position/genes do
not coincide with reported positions/gene, select the \"No valid match\"
option. This will initiate a non-exact search.""",
            '5': """For each unmatched sites in the previous form, a soft search was performed and
possible matches are being displayed.  If none of possible matches are good
enough, you can prefer not matching by selecting the last option for a site.
Inexact matches for sites without valid matches are listed here, sorted by
affinity to the TF-binding motif.  If the matched position/genes do not coincide
with reported positions/gene, select the \"No valid match\" option.""",
            '6': """For each transcription factor binding site positioned in the genome, nearby
genes are being displayed. If the curated publication contains any experimental
support for any TF-gene regulation, mark those genes. All other genes will be
saved as 'inferred' regulation. If the publication was marked as it doesn't
contain expression data, checkboxes should be disabled. Nearby genes are
displayed for identified sites. Check all genes for which TF-site mediated
regulation is reported in the manuscript. Skip this step if manuscript does not
report gene expression.""",
            '7': "This is the last form before submission of the curation. " \
                 "Fill all required fields."
        }

        context["form_title"] = titles[self.steps.current]
        context["form_description"] = descriptions[self.steps.current]
        return context

    def get_form(self, step=None, data=None, files=None):
        """Construct the form for a given step. The method is overriden to
        add/update arguments to the form instance."""
        form = super(CurationWizard, self).get_form(step, data, files)
        if step == None:
            step = self.steps.current

        handlers = {'0': publication_get_form, # functions
                    '1': genome_get_form,
                    '2': techniques_get_form,
                    '3': site_report_get_form,
                    '4': site_exact_match_get_form,
                    '5': site_soft_match_get_form,
                    '6': site_regulation_get_form,
                    '7': curation_review_get_form,
                   }
        
        return handlers[step](wiz=self, form=form)

    def process_step(self, form):
        """Process data after each step before it gets stored"""
        handlers = {'0': publication_process, # functions
                    '1': genome_process,
                    '2': techniques_process,
                    '3': site_report_process,
                    '4': site_exact_match_process,
                    '5': site_soft_match_process,
                    '6': site_regulation_process,
                    '7': curation_review_process,
                   }
        
        handlers[self.steps.current](wiz=self, form=form)

        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        """Last step in curation process, this method is called after all
        forms. Insert all data into the database."""
        # get data from forms
        publication = publication_done(self, form_list[0], **kwargs)
        genome_cd = genome_done(self, form_list[1], **kwargs)
        techniques_cd = techniques_done(self, form_list[2], **kwargs)
        curation_review_cd = curation_review_done(self, form_list[7], **kwargs)

        # find curator
        curator = models.Curator.objects.get(user=self.request.user)
        # create curation object
        curation = models.Curation(publication=publication,
                                   TF_species=genome_cd['TF_species'],
                                   site_species=genome_cd['site_species'],
                                   TF=genome_cd['TF'],
                                   TF_function=genome_cd['TF_function'],
                                   TF_type=genome_cd['TF_type'],
                                   TF_instance=genome_cd['TF_instance'],
                                   experimental_process=techniques_cd['experimental_process'],
                                   forms_complex=techniques_cd['forms_complex'],
                                   complex_notes=techniques_cd['complex_notes'],
                                   requires_revision=curation_review_cd['requires_revision'],
                                   notes=curation_review_cd['notes'],
                                   confidence=curation_review_cd['confidence'],
                                   NCBI_submission_ready=curation_review_cd['NCBI_submission_ready'],
                                   curator=curator)
        curation.save()

        # add techniques
        for t in techniques_cd['techniques']:
            curation.experimental_techniques.add(t)

        regulations = site_regulation_done(self, form_list[6], **kwargs)
        # create matched site instances & regulations
        site_match_done(self, curation, regulations)
        # create not annotated sites objects
        not_matched_sites_done(self, curation)

        # mark paper as complete if so
        if curation_review_cd['paper_complete']:
            publication.curation_complete = True
            publication.save()
            
        # if this curation finalizes the curation of the paper, mark paper as so.
        if sutils.sin(self.request.session, 'old_curation'):
            # delete existing one
            sutils.sget(self.request.session, 'old_curation').delete()
        
        return HttpResponseRedirect(reverse(views.success))
    
    
# curation handler

# for form definitions, go curationform.py
@login_required
def curation(request):
    # clear session data from previous forms
    
    # If user selects the old curation and then go back, the session will have the
    # old_curation key in table, and it will cause trouble.
    if sutils.sin(request.session, 'old_curation'):
        sutils.sdel(request.session, 'old_curation')
    # TODO make custom session objects be form wizard based
    
    view = CurationWizard.as_view([PublicationForm,
                                   GenomeForm,
                                   TechniquesForm,
                                   SiteReportForm,
                                   SiteExactMatchForm,
                                   SiteSoftMatchForm,
                                   SiteRegulationForm,
                                   CurationReviewForm])
    return view(request)
                                   
