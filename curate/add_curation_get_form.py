"""This file contains get_form functions for curation-wizard steps. At each step
of the curation, get_form function is called to construct the form for the
particular step.

get_form documentation:

https://docs.djangoproject.com/en/1.6/ref/contrib/formtools/form-wizard/#django.contrib.formtools.wizard.views.WizardView.get_form

"""

import models
import templatetags.print_pub
import session_utils
from django import forms
from django.utils.safestring import mark_safe
from django.contrib import messages


def publication_get_form(wiz, form):
    """Construct the form for publication selection step."""
    user = wiz.request.user
    curator = models.Curator.objects.get(user=user)
    # select papers assigned to user and not complete
    assigned_pubs = models.Publication.objects.filter(assigned_to=curator,
                                                      curation_complete=False)
    # put them in form choices, populate form field
    choices = [(p.publication_id, mark_safe("%s" % (templatetags.print_pub.print_pub(p))))
               for p in assigned_pubs]
    form.fields["pub"].choices = choices
    return form


def genome_get_form(wiz, form):
    """Construct the form for genome and TF selection step."""
    c = session_utils.get(wiz.request.session, "previously_curated_paper")
    # If selected publication is the one most recently curated, the related curation
    # should be in object c. Otherwise, c = None.  If so, populate "Genome and TF
    # information" form fields from the previously submitted curation to make things
    # easier for curator.
    if c:
        form.initial["TF"] = c.TF
        form.initial["TF_type"] = c.TF_type
        if c.site_instances.all():
            form.initial["genome_accession"] = c.site_instances.all()[0].genome.genome_accession
            if len(c.site_instances.all()) == 2:
                form.initial['genome_accession_1'] = c.site_instances.all()[1].genome.genome_accession
            if len(c.site_instances.all()) == 3:
                form.initial['genome_accession_1'] = c.site_instances.all()[2].genome.genome_accession
            
        # Enter TF accession numbers
        form.initial["TF_accession"] = c.TF_instances.all()[0].protein_accession
        if len(c.TF_instances.all()) == 2:
            form.initial["TF_accession_1"] = c.TF_instances.all()[1].protein_accession
        if len(c.TF_instances.all()) == 3:
            form.initial["TF_accession_2"] = c.TF_instances.all()[2].protein_accession
        
        form.initial["TF_species"] = c.TF_species
        form.initial["site_species"] = c.site_species

        msg = """
        <h4>Warning!</h4> It seems that the paper you selected is previously
        curated. For convenience, fields in this form are automatically filled based on
        the previous curation of the paper. They may differ in this curation, so it is
        best to check that they are correct before proceeding to the next step."""
        messages.warning(wiz.request, mark_safe(msg))

    # In addition populate two fields on whether the manuscript contains
    # experimental data and promoter information
    pid = session_utils.get(wiz.request.session, 'publication')
    pub = models.Publication.objects.get(publication_id=pid)
    form.fields["contains_promoter_data"].initial = pub.contains_promoter_data
    form.fields["contains_expression_data"].initial = pub.contains_expression_data

    return form

def techniques_get_form(wiz, form):
    """Construct the form for experiemental techniques step."""
    c = session_utils.get(wiz.request.session, 'previously_curated_paper')
    # if selected paper is previously curated, prepopulate experimental techniques
    if c:
        # get all techniques used in this curation
        cur_site_insts = models.Curation_SiteInstance.objects.filter(curation=c)
        techniques = list(set(t.technique_id for csi in cur_site_insts
                              for t in csi.experimental_techniques.all()))
        form.fields['techniques'].initial = techniques
        
        form.fields['experimental_process'].initial = c.experimental_process
        try:
            external_db = models.Curation_ExternalDatabase.objects.get(curation=c)
        except models.Curation_ExternalDatabase.DoesNotExist:
            external_db = None
        form.fields['external_db_type'].initial = external_db.external_database.ext_database_id if external_db else None
        form.fields['external_db_accession'].initial = external_db.accession_number if external_db else ""
        form.fields['forms_complex'].initial = c.forms_complex
        form.fields['complex_notes'].initial = c.complex_notes
    return form

def site_entry_get_form(wiz, form):
    """Construct the form for site entry step."""
    c = session_utils.get(wiz.request.session, 'previously_curated_paper')
    # if paper is previously curated, prepopulate fields
    if c:
        # pick any curation_site_instance object for this curation
        try:
            curation_site_instance = models.Curation_SiteInstance.objects.filter(curation=c).all()[:1].get()
            #form.fields['is_motif_associated'].initial = curation_site_instance.is_motif_associated
            form.fields['site_type'].initial = curation_site_instance.site_type
        except:
            pass
        # Delete session data, if user change any field and then come back,
        # Store users last entered data, instead of populated data.
        session_utils.put(wiz.request.session, "previously_curated_paper", None)
        
    return form

def site_exact_match_get_form(wiz, form):
    """Show the list of sites and their exact matches in this form.  Each
    reported site and its exact match results are represented as a field. The
    form is generated dynamically as the number of exact matches is not
    fixed."""
    sites = session_utils.get(wiz.request.session, "sites")
    for i,site in enumerate(sites):
        label = site.seq
        choices = site.populate_match_choices(add_no_valid_opt=True)
        form.fields[i] = forms.ChoiceField(label=label, choices=choices,
                                            widget=forms.RadioSelect())
        form.fields[i].initial = str(choices[0][0])
            
    return form

def site_soft_match_get_form(wiz, form):
    """For sites that are not matched in the previous step (exact match step),
    show soft-search results."""
    sites = session_utils.get(wiz.request.session, "sites")
    for i,site in enumerate(sites):
        if not site.is_matched():
            label = site.seq
            choices = site.populate_match_choices(add_no_valid_opt=True)
            form.fields[i] = forms.ChoiceField(label=label, choices=choices,
                                               widget=forms.RadioSelect())
            form.fields[i].initial = str(choices[0][0])
    return form

def site_annotation_get_form(wiz, form):
    """Annotation step for site instances."""
    sites = session_utils.get(wiz.request.session, 'sites')
    techniques = session_utils.get(wiz.request.session, 'techniques')
    for i,site in enumerate(sites):
        if site.is_matched():
            # create dummy field for the label
            form.fields['%d_site' % i] = forms.BooleanField(label=site.get_match().pprint(), required=False)
            # create quantitative value field
            if session_utils.get(wiz.request.session, 'has_quantitative_data'):
                form.fields['%d_qval' % i] = forms.FloatField(label='Q-val',
                                                              required=False,
                                                              initial=site.qval)
            # create TF function selection field
            form.fields['%d_TF_function' % i] = forms.ChoiceField(label='TF function',
                                                                  choices=models.Curation.TF_FUNCTION)
            # create techniques fields (one checkbox for each technique)
            for j,t in enumerate(techniques):
                form.fields['%d_technique_%d' % (i,j)] = forms.BooleanField(label="", required=False)
            
    return form

def gene_regulation_get_form(wiz, form):
    """For each site, show nearby genes for regulation input. For each nearby
    gene close to a site, curator manually checks if site regulates the gene or
    not (i.e. if there is any experimental support for regulation)."""
    # If publication is marked as not having expression data, disable the form
    # after building it.
    pid = session_utils.get(wiz.request.session, 'publication')
    publication = models.Publication.objects.get(pk=pid)
    # Get all site matches
    sites = session_utils.get(wiz.request.session, 'sites')
    for i,site in enumerate(sites):
        if site.is_matched():
            choices = [(g.gene_id, "%s (%s): %s" % (g.locus_tag, g.name, g.description))
                       for g in site.get_match().nearby_genes]
            form.fields[i] = forms.MultipleChoiceField(label=site.seq, choices=choices,
                                                       required=False,
                                                       widget=forms.CheckboxSelectMultiple(),
                                                       help_text = site.get_match().match_diagram)
            if not publication.contains_expression_data:
                form.fields[i].widget.attrs['disabled'] = 'disabled'
                
    return form

def curation_review_get_form(wiz, form):
    return form
