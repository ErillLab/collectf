"""

Handler class for curation process.

From Django docs:

Here's the basic workflow for how a user would use a wizard:

1) The user visits the first page of the wizard, fills in the form and submits it.
n2) The server validates the data. If it's invalid, the form is displayed again,
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

from django.utils.safestring import mark_safe
from django.contrib import messages
from baseapp.templatetags import gene_diagram
from baseapp.templatetags import publication_tags
from baseapp.templatetags import pretty_print

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
    choices = [(p.publication_id, mark_safe("[%s] %s" % (p.publication_id, publication_tags.as_publication(p))))
               for p in not_completed_pubs]
    form.fields["pub"].choices = choices
    return form
    
def genome_get_form(wiz, form):
    """Genome, TF, TF_family, tf instance, .. selection form"""

    c = sutils.sget(wiz.request.session, "previously_curated_paper")
    # If selected publication is the one most recently curated, the related curation
    # should be in object c. Otherwise, c = None.  If so, populate "Genome and TF
    # information" form fields from the previously submitted curation to make things
    # easier for curator.
    if c:
        form.initial["TF"] = c.TF
        form.initial["TF_type"] = c.TF_type
        form.initial["TF_function"] = c.TF_function
        if c.site_instances.all():
            form.initial["genome_accession"] = c.site_instances.all()[0].genome.genome_accession
        form.initial["TF_accession"] = c.TF_instance.protein_accession
        
        form.initial["TF_species"] = c.TF_species
        form.initial["site_species"] = c.site_species

        msg = """
        <h4>Warning!</h4> It seems that the paper you selected is previously
        curated. For convenience, fields in this form are automatically filled based on
        the previous curation of the paper. They may differ in this curation, so it is
        best to check that they are correct before proceeding to the next step."""
        
        messages.warning(wiz.request, mark_safe(msg))

    return form

def techniques_get_form(wiz, form):
    # populate the following fields
    # - contains expression data
    # - contains promoter data
    pid = sutils.sget(wiz.request.session, 'publication')
    pub = models.Publication.objects.get(publication_id=pid)
    form.fields["contains_promoter_data"].initial = pub.contains_promoter_data
    form.fields["contains_expression_data"].initial = pub.contains_expression_data
    c = sutils.sget(wiz.request.session, 'previously_curated_paper')
    # if selected paper is previously curated, prepopulate experimental techniques
    if c:
        print c.curation_id
        form.fields['techniques'].initial = [str(t.technique_id) for t in c.experimental_techniques.all()]
        form.fields['experimental_process'].initial = c.experimental_process
        try:
            external_db = models.Curation_ExternalDatabase.objects.get(curation=c)
        except models.Curation_ExternalDatabase.DoesNotExist:
            external_db = None
        form.fields['external_db_type'].initial = external_db.external_database.ext_database_id if external_db else None
        form.fields['external_db_accession'].initial = external_db.accession_number if external_db else ""
        form.fields['forms_complex'].initial = c.forms_complex
        form.fields['complex_notes'].initial = c.complex_notes

        # delete session data, if user change any field and then come back,
        # store users last entered data, instead of populated data.
        sutils.sput(wiz.request.session, "previously_curated_paper", None)

    return form

def site_report_get_form(wiz, form):
    #msg = "hello site report form"
    #messages.info(wiz.request, msg)
    return form

# helper function for site_exact_match_form and site_soft_match_forms
def populate_match_choices(site, matches, is_exact, add_no_valid_opt=False):
    # for a given site and its all matches, populate Django field choice
    # exact_match: True if populating exact match results, if false, soft search
    choices = []
    for mid, match in matches.items():
        choices.append((mid, pretty_print.print_site_match(site, match, is_exact)))
    if add_no_valid_opt:
        last_choice_msg = "No valid match"
        choices.append((None, last_choice_msg))
    return choices
    
def site_exact_match_get_form(wiz, form):
    """Show list of sites and their exact matches."""
    # Each reported site and its exact match results are represented as field.
    # Dynamically add fields to site search results step.
    sites = sutils.sget(wiz.request.session, 'sites')
    site_match_choices = sutils.sget(wiz.request.session, 'site_match_choices')

    is_coordinate = sutils.sget(wiz.request.session, 'is_coordinate')

    # If sites are reported as coordinates, this form is read only
    add_no_valid_opt = not is_coordinate
    for sid, matches in site_match_choices.items(): # for all matches belong to a site
        label = pretty_print.site2label(sid,sites[sid])
        choices = populate_match_choices(sites[sid],
                                         matches,
                                         is_exact=True,
                                         add_no_valid_opt=add_no_valid_opt)
        # make the form field
        form.fields[sid] = forms.ChoiceField(label=label,
                                             choices=choices,
                                             widget=forms.RadioSelect())
        

        if is_coordinate: # make the first option selected
            #form.fields[sid].widget.attrs['disabled']='disabled'
            form.fields[sid].initial = choices[0][0]
            
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
        label = pretty_print.site2label(sid, sites[sid])
        choices = populate_match_choices(sites[sid], matches, is_exact=False, add_no_valid_opt=True)
        # make the form field
        form.fields[sid] = forms.ChoiceField(label=label,
                                             choices=choices,
                                             widget=forms.RadioSelect())
    return form

def site_quantitative_data_get_form(wiz, form):
    sites = sutils.sget(wiz.request.session, 'sites')
    exact_site_matches = sutils.sget(wiz.request.session, 'exact_site_matches')
    soft_site_matches = sutils.sget(wiz.request.session, 'soft_site_matches')

    all_site_matches = []
    all_site_matches += exact_site_matches.items() if exact_site_matches else []
    all_site_matches += soft_site_matches.items() if soft_site_matches else []
    all_site_matches = dict(all_site_matches)
    
    site_quantitative_data = sutils.sget(wiz.request.session, 'site_quantitative_data')
    assert set(sites) == set(site_quantitative_data) # make sure keys are same
    for sid in all_site_matches:
        form.fields[sid] = forms.FloatField(
            label=pretty_print.site2label(sid,sites[sid]),
            help_text=pretty_print.print_site_match(sites[sid], all_site_matches[sid], is_exact=True),
            initial=site_quantitative_data[sid])

    print 'form_ready'
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

    all_site_matches = []
    all_site_matches += exact_site_matches.items() if exact_site_matches else []
    all_site_matches += soft_site_matches.items() if soft_site_matches else []
    
    for sid,match in all_site_matches:
        choices = []  # list of genes to the site match
        for g in match.nearby_genes:
            choices.append((g.gene_id, '%s (%s): %s' % (g.locus_tag, g.name, g.description)))

        label = pretty_print.match2label(sid, match)

        form.fields[sid] = forms.MultipleChoiceField(label=label,
                                                     choices=choices,
                                                     required=False,
                                                     widget=forms.CheckboxSelectMultiple(),
                                                     help_text=gene_diagram.site_match_diagram(match))
        
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

    if form.cleaned_data["no_data"]:
        # mark paper as having no data
        paper = models.Publication.objects.get(publication_id=pubid)
        note = " \nPaper has no TF-binding site data."
        paper.submission_notes += note
        paper.curation_complete = True
        paper.save()
        sutils.sput(wiz.request.session, "paper_contains_no_data", True)
        return


    # if paper is previously curated, populate genome and TF information form
    # search db if there is any curation object belonging to this publication
    p = models.Publication.objects.get(publication_id=pubid)
    # check if the publication is previously curated
    cs = models.Curation.objects.filter(publication=p)
    if cs.count() >= 1:
        sutils.sput(wiz.request.session, "previously_curated_paper", cs[0])
    else:
        sutils.sput(wiz.request.session, "previously_curated_paper", None)

def genome_process(wiz, form):
    genome_accession = form.cleaned_data['genome_accession']
    # in form validation genome is searched in db, and if not found,
    # it is inserted into db. So, at this point, it is guaranteed that
    # genome with id <genome_accession> should be in db.
    genome = models.Genome.objects.get(genome_accession=genome_accession)
    # store genome in session data
    sutils.sput(wiz.request.session, 'genome', genome)
    # store site species in session data
    sutils.sput(wiz.request.session, 'site_species', form.cleaned_data['site_species'])
    

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
    """For site report form, the user may enter either site sequences or
    coordinates. In addition to that, quantitative data associated with sites may be
    given as well. The case that makes everything complicated (and makes me writing
    this comment at 1 am in the morning) is that user may prefer to enter those
    quantitative values in another form step (after the site match step)."""
    def process_helper_sites_only():
        # process data in case of only site sequence data
        # Get reported sites -- [(sid, site), ..]
        sites = sitesearch.parse_site_input(form.cleaned_data['sites'].strip())
        # Find exact matches
        site_match_choices = sitesearch.match_all_exact(genome, genes, sites)
        # Put sites and site matches into session data
        sutils.sput(wiz.request.session, 'site_match_choices', site_match_choices)
        sutils.sput(wiz.request.session, 'sites', sites)
        # no quantitative data for this case
        sutils.sput(wiz.request.session, 'site_quantitative_data', {})

    def process_helper_sites_with_quantitative_data():
        # process data of site sequences with quantitative values
        # parse sites and quantitative values
        lines = re.split('[\r\n]+', form.cleaned_data['sites'].strip())
        sites = [line.split()[0] for line in lines]
        quantitative_values = [line.split()[1] for line in lines]
        assert len(sites) == len(quantitative_values)
        # Get reported sites
        sites = sitesearch.parse_site_input('\n'.join(sites))
        site_quantitative_data = dict(enumerate(map(float, quantitative_values)))
        # Find exact matches
        site_match_choices = sitesearch.match_all_exact(genome, genes, sites)
        # Put sites and site matches into session data
        sutils.sput(wiz.request.session, 'site_match_choices', site_match_choices)
        sutils.sput(wiz.request.session, 'sites', sites)
        sutils.sput(wiz.request.session, 'site_quantitative_data', site_quantitative_data)
        print sites
        print site_quantitative_data
        
    def process_helper_coordinates(with_quantitative=False):
        # process data in case of only coordinate data is given.
        coordinates = form.cleaned_data['sites'].strip()
        coordinates = [re.split('[\t ]+', line) for line in re.split('[\r\n]+', coordinates)]
        res = sitesearch.match_all_exact_coordinates_only(genome, genes, coordinates)
        sites, site_match_choices, site_quantitative_data = res
        # This assertion should stay here (site_quantitative_data = {} -> request.session below)
        if not with_quantitative: # no quantitative data should be present
            assert all(not val for val in site_quantitative_data.values()), "site_quantitative_data not null"
            site_quantitative_data = {}
        # Put sites and site matches into session data
        sutils.sput(wiz.request.session, 'site_match_choices', site_match_choices)
        sutils.sput(wiz.request.session, 'sites', sites)
        sutils.sput(wiz.request.session, 'site_quantitative_data', site_quantitative_data)


    def process_helper_chip_assoc():
        # The list of motif-associated sites are read from -sites- field.  The list
        # of ChIP sequences (not-motif-associated-sites), with peak-intensity-values
        # are read from ChIP-extra-field and peak-intensity-values are associated
        # with motif-associated-sites.

        print 'process_helper_chip_assoc'
        # make sure quantitative values exist.
        assert form.cleaned_data['has_quantitative_data'], "has_quantitative_data error"
        # First, check whether <sites> field in sequence or coordinate format.
        if form.cleaned_data['is_coordinate']:
            process_helper_coordinates(with_quantitative=False)
        else:
            process_helper_sites_only()

        # process chip-extra-field
        coords = form.cleaned_data['chip_data_extra_field'].strip()
        coords = [re.split('[\t ]+', line) for line in re.split('[\r\n]+', coords)]
        # this is just to get sequences for coordinates in Chip_data_extra_field
        chip_peaks,_,quantitative_vals = sitesearch.match_all_exact_coordinates_only(genome, genes, coords)
        
        assert not sutils.sget(wiz.request.session, 'site_quantitative_data')
        # Now, associate quantitative_data (peak_intensity values) with site instances
        sites = sutils.sget(wiz.request.session, 'sites')
        site_quantitative_data = {}
        for site_id, site in sites.items():
            # search that site in ChIP extra field
            for peak_id, peak_seq in chip_peaks.items():
                if site in peak_seq or sitesearch.reverse_complement(site) in peak_seq:
                    site_quantitative_data[site_id] = quantitative_vals[peak_id]
                    break
            else:
                site_quantitative_data[site_id] = None
                
        assert len(site_quantitative_data) == len(sites)
        sutils.sput(wiz.request.session, 'site_quantitative_data', site_quantitative_data)
        sutils.sput(wiz.request.session, 'soft_site_matches', {})
        sutils.sput(wiz.request.session, 'not_matched_sites', [])
        sutils.sput(wiz.request.session, 'soft_site_match_choices', {})
    
    def site_report_process_helper_0():
        #p NOT is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_only()

    def site_report_process_helper_1():
        # NOT is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # is_coordinate
        process_helper_coordinates(with_quantitative=False)

    def site_report_process_helper_2():
        # NOT is_motif_associated
        # NOT is_chip_data
        # has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_with_quantitative_data()

    def site_report_process_helper_3():
        # NOT is_motif_associated
        # NOT is_chip_data
        #  has_quantitative_data
        #  is_coordinate
        process_helper_coordinates(with_quantitative=True)
        
    def site_report_process_helper_4():
        # NOT is_motif_associated
        #  is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_only()

    def site_report_process_helper_5():
        # NOT is_motif_associated
        #  is_chip_data
        # NOT has_quantitative_data
        #  is_coordinate
        process_helper_coordinates(with_quantitative=False)

    def site_report_process_helper_6():
        # NOT is_motif_associated
        #  is_chip_data
        #  has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_with_quantitative_data()

    def site_report_process_helper_7():
        # NOT is_motif_associated
        #  is_chip_data
        #  has_quantitative_data
        #  is_coordinate
        process_helper_coordinates(with_quantitative=True)

    def site_report_process_helper_8():
        # is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_only()
        
    def site_report_process_helper_9():
        #  is_motif_associated
        # NOT is_chip_data
        # NOT has_quantitative_data
        #  is_coordinate
        process_helper_coordinates(with_quantitative=False)

    def site_report_process_helper_10():
        #  is_motif_associated
        # NOT is_chip_data
        #  has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_with_quantitative_data()

    def site_report_process_helper_11():
        #  is_motif_associated
        # NOT is_chip_data
        #  has_quantitative_data
        #  is_coordinate
        process_helper_coordinates(with_quantitative=True)

    def site_report_process_helper_12():
        #  is_motif_associated
        #  is_chip_data
        # NOT has_quantitative_data
        # NOT is_coordinate
        process_helper_sites_only()

    def site_report_process_helper_13():
        #  is_motif_associated
        #  is_chip_data
        # NOT has_quantitative_data
        #  is_coordinate
        process_helper_coordinates(with_quantitative=False)
        
    def site_report_process_helper_14():
        #  is_motif_associated
        #  is_chip_data
        #  has_quantitative_data
        # NOT is_coordinate
        # SPECIAL CASE
        process_helper_chip_assoc()

    def site_report_process_helper_15():
        #  is_motif_associated
        #  is_chip_data
        #  has_quantitative_data
        #  is_coordinate
        process_helper_chip_assoc()

    print 'process_begin'

    sutils.sput(wiz.request.session, 'sites', {})
    sutils.sput(wiz.request.session, 'exact_site_matches', {})
    sutils.sput(wiz.request.session, 'soft_site_matches', {})
    sutils.sput(wiz.request.session, 'not_matched_sites', [])
    sutils.sput(wiz.request.session, 'site_quantitative_data', {})
    sutils.sput(wiz.request.session, 'soft_site_match_choices', {})

    is_motif_associated = form.cleaned_data.get('is_motif_associated')
    is_chip_data = form.cleaned_data.get('is_chip_data')
    is_coordinate = form.cleaned_data.get('is_coordinate')
    has_quantitative_data = form.cleaned_data.get('has_quantitative_data')

    # read genome from session data
    genome = sutils.sget(wiz.request.session, 'genome')
    genes = models.Gene.objects.filter(genome=genome).order_by('start')

    call_func_id = (('1' if is_motif_associated else '0') +
                    ('1' if is_chip_data else '0') +
                    ('1' if has_quantitative_data else '0') +
                    ('1' if is_coordinate else '0'))
    
    call_func_str = 'site_report_process_helper_%d()' % int(call_func_id,2)
    print 'process_func', call_func_str
    eval(call_func_str)



    # store booleans
    sutils.sput(wiz.request.session, 'is_motif_associated', is_motif_associated)
    sutils.sput(wiz.request.session, 'is_chip_data', is_chip_data)
    sutils.sput(wiz.request.session, 'is_coordinate', is_coordinate)
    sutils.sput(wiz.request.session, 'has_quantitative_data', has_quantitative_data)

    if form.cleaned_data.get('has_quantitative_data'):
        print form.cleaned_data.get('quantitative_data_format')
        sutils.sput(wiz.request.session, 'quantitative_data_format', form.cleaned_data.get('quantitative_data_format'))
    else:
        sutils.sput(wiz.request.session, 'site_quantitative_data', {})
        sutils.sput(wiz.request.session, 'quantitative_data_format', None)

    if is_chip_data:
        sutils.sput(wiz.request.session, 'assay_conditions', form.cleaned_data.get('assay_conditions'))
        sutils.sput(wiz.request.session, 'chip_method_notes', form.cleaned_data.get('chip_method_notes'))
    else:
        sutils.sput(wiz.request.session, 'assay_conditions', None)
        sutils.sput(wiz.request.session, 'chip_method_notes', None)
    
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


def site_quantitative_data_process(wiz, form):
    quantitative_vals = sutils.sget(wiz.request.session, 'site_quantitative_data')
    for id,val in form.cleaned_data.items():
        quantitative_vals[id] = val
    sutils.sput(wiz.request.session, 'site_quantitative_data', quantitative_vals)
        
def site_regulation_process(wiz, form):
    pass

def curation_review_process(wiz, form):
    pass

# CurationWizard done step functions

# When the last step of the form is submitted, gather all data from each form
# step and insert into database. Following _done_ functions are implemented for
# each step and they are called by CurationWizard.done when user submit curation
# form.
def publication_done(wiz):
    pid = sutils.sget(wiz.request.session, 'publication')
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
    cd['external_db_type'] = form.cleaned_data['external_db_type']
    cd['external_db_accession'] = form.cleaned_data['external_db_accession']
    return cd

def site_report_done(wiz, form, **kwargs):
    pass

def site_exact_match_done(wiz, form, **kwargs):
    """Get data from exact site match form"""
    pass

def site_soft_match_done(wiz, form, **kwargs):
    """Get data from soft site match form"""
    pass

def site_quantitative_data_done(wiz, form, **kwargs):
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
    quantitative_vals = sutils.sget(wiz.request.session, 'site_quantitative_data')
    genome = sutils.sget(wiz.request.session, 'genome')  # get genome
    # combine exact and soft site matches
    all_site_matches = []
    all_site_matches += exact_site_matches.items() if exact_site_matches else []
    all_site_matches += soft_site_matches.items() if soft_site_matches else []
    all_site_matches = dict(all_site_matches)

    is_motif_associated = sutils.sget(wiz.request.session, 'is_motif_associated')
    

    for sid, match in all_site_matches.items():
        # create SiteInstance object (or get if available)
        si, created  = models.SiteInstance.objects.get_or_create(
            genome=genome,
            start=match.match.start,
            end=match.match.end,
            strand=match.match.strand,
            seq=match.match.seq)

                
        # create curation_siteinstance object (through relation)
        cs = models.Curation_SiteInstance(curation=curation,
                                          site_instance=si,
                                          annotated_seq=sites[sid],
                                          is_motif_associated=is_motif_associated,
                                          #chip_info=chip_data,
                                          #quantitative_data_format=quantitative_data_format,
                                          quantitative_value = quantitative_vals.get(sid, None))
        cs.save()

        # for site instance, add regulation information
        regulation_done(wiz,
                        regulations=regulations[sid],
                        match=match,
                        curation_site_instance=cs)


        
def not_matched_sites_done(wiz, curation):
    """Create not matched site instances."""
    sites = sutils.sget(wiz.request.session, 'sites')
    not_matched_sites = sutils.sget(wiz.request.session, 'not_matched_sites')
    if not_matched_sites:
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
            '6': "Site-Quantitative Value Association",
            '7': "Gene regulation (experimental support)",
            '8': "Curation information"
        }
        descriptions = {
            '0': "Please choose a publication to curate.",
            '1': """This step collects information on the transcription factor (TF), the specific
            strains reported in the manuscript and the NCBI GenBank sequences that reported
            sites and TF will be mapped onto.""",
            '2': """Select the experimental techniques and the describe the basic experimental
            procedure used to verify binding/expression of the sites reported in this
            curation.""",
            '3': '', # filled dynamically
            '4': """For each reported site, all exact matches in the chosen genome are listed. If a
            reported site does not have any exact matches, or the matched position/genes do
            not coincide with reported positions/gene, select the \"No valid match\"
            option. This will initiate a non-exact search.""",
            '5': """Inexact matches for sites without valid matches are listed here, sorted by
            affinity to the TF-binding motif.  If the matched position/genes do not coincide
            with reported positions/gene, select the \"No valid match\" option.""",
            '6': "description goes here",
            '7': """Nearby genes are
            displayed for identified sites. Check all genes for which TF-site mediated
            regulation is reported in the manuscript. Skip this step if manuscript does not
            report gene expression.""",
            '8': """This step finalizes the curation. Fill all required fields."""
        }
        context["form_title"] = titles[self.steps.current]
        context["form_description"] = descriptions[self.steps.current]

        # add some extra content for some steps

        if self.steps.current == '4':
            sites = sutils.sget(self.request.session, 'sites')
            is_motif_associated = sutils.sget(self.request.session, 'is_motif_associated')
            if is_motif_associated and len(sites) > 1:
                context.update({'weblogo_img': bioutils.weblogo_uri(sites.values())})
        
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
                    '6': site_quantitative_data_get_form,
                    '7': site_regulation_get_form,
                    '8': curation_review_get_form,
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
                    '6': site_quantitative_data_process,
                    '7': site_regulation_process,
                    '8': curation_review_process,
                   }
        
        handlers[self.steps.current](wiz=self, form=form)

        return self.get_form_step_data(form)

    def render(self, form=None, **kwargs):
        form = form or self.get_form()
        context = self.get_context_data(form=form, **kwargs)

        # the only reason to override this method is to check "no data" button
        # in publication selection form step. If the curator selects a
        # publication and checks the button "This paper contains no data.", the
        # proper action is to terminate curation process as ther is no data to
        # be curated. The problem is to process that information and redirect to
        # the home page. This is achieved by overloading render method. Before
        # this function is called, publication object is modified as having no
        # data. Afterwards, this function is called and it redirects to the
        # homepage with the message about the action that was performed.

        if (sutils.sin(self.request.session, "paper_contains_no_data") and
            sutils.sget(self.request.session, "paper_contains_no_data")):
            msg = "The publication was marked as having no data."
            messages.info(self.request, msg)
            # clear session data
            sutils.clear(self.request.session)
            return HttpResponseRedirect(reverse(views.home))
        
        return self.render_to_response(context)

    def done(self, form_list, **kwargs):
        """Last step in curation process, this method is called after all
        forms. Insert all data into the database."""
        # get data from forms
        publication = publication_done(self)
        
        # quick and dirty
        # If curation view is used for revising a curation, the first step
        # (publication selection) is hidden in formwizard, and in form_list
        # list. However, in that case, the first element (i.e. form_list[0]) is
        # not about publication anymore. Therefore, the quick and dirty solution
        # is to insert a "fake" element to the form_list list
        if sutils.sin(self.request.session, 'old_curation'):
             # this is a revision for an existing curation
             form_list.insert(0, None)

        head = lambda x: x[0]

        genome_form = head([f for f in form_list if type(f) == GenomeForm])
        genome_cd = genome_done(self, genome_form, **kwargs)

        techniques_form = head([f for f in form_list if type(f) == TechniquesForm])
        techniques_cd = techniques_done(self, techniques_form, **kwargs)

        curation_review_form = head([f for f in form_list if type(f) == CurationReviewForm])
        curation_review_cd = curation_review_done(self, curation_review_form, **kwargs)

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

        is_chip_data = sutils.sget(self.request.session, 'is_chip_data')
        has_quantitative_data = sutils.sget(self.request.session, 'has_quantitative_data')
        chip_data = None
        quantitative_data_format = None
        
        if is_chip_data:
            assay_conditions = sutils.sget(self.request.session, 'assay_conditions')
            chip_method_notes = sutils.sget(self.request.session, 'chip_method_notes')
            chip_data = models.ChipInfo(assay_conditions=assay_conditions,
                                        method_notes=chip_method_notes)
            chip_data.save()

        if has_quantitative_data:
            quantitative_data_format = sutils.sget(self.request.session, 'quantitative_data_format')
            print 'quant format', quantitative_data_format

        curation.chip_info = chip_data
        curation.quantitative_data_format = quantitative_data_format
        curation.save()

        # add techniques
        for t in techniques_cd['techniques']:
            curation.experimental_techniques.add(t)

        # add external db references
        if techniques_cd['external_db_type'] != "None" and techniques_cd['external_db_accession']:
            external_db_type = models.ExternalDatabase.objects.get(ext_database_id=techniques_cd['external_db_type'])
            curation_ext_ref = models.Curation_ExternalDatabase(curation = curation,
                                                                external_database=external_db_type,
                                                                accession_number=techniques_cd['external_db_accession'])
            curation_ext_ref.save()

        curation.save()

        site_regulation_form = head([f for f in form_list if type(f) == SiteRegulationForm])
        regulations = site_regulation_done(self, site_regulation_form, **kwargs)
        # create matched site instances & regulations
        site_match_done(self, curation, regulations)
        # create not annotated sites objects
        not_matched_sites_done(self, curation)

        # mark paper as complete if so
        if curation_review_cd['paper_complete']:
            publication.curation_complete = True
            publication.save()

        else:
            publication.curation_complete = False
            publication.save()
            
        # if this curation finalizes the curation of the paper, mark paper as so.
        if sutils.sin(self.request.session, 'old_curation'):
            # delete existing one
            sutils.sget(self.request.session, 'old_curation').delete()

        # clear session
        sutils.clear(self.request.session)
        
        messages.success(self.request, "Curation was successfully submitted.")
        return HttpResponseRedirect(reverse(views.home))




# curation handler

# for form definitions, go curationform.py
@login_required
def curation(request):

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
                                   SiteQuantitativeDataForm,
                                   SiteRegulationForm,
                                   CurationReviewForm],
                                  condition_dict = {'4': exact_site_match_form_condition,
                                                    '5': inexact_site_match_form_condition,
                                                    '6': site_quantitative_data_form_condition})

    return view(request)

def exact_site_match_form_condition(wizard):
    return True

def inexact_site_match_form_condition(wizard):
    return sutils.sget(wizard.request.session, 'soft_site_match_choices')

def site_quantitative_data_form_condition(wizard):
    return sutils.sget(wizard.request.session, 'has_quantitative_data')
