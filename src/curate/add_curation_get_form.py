"""This module contains get_form functions for curation-wizard steps. At each step
of the curation, get_form function is called to construct the form for the
particular step.
"""

from django import forms
from django.utils.safestring import mark_safe
from django.contrib import messages

from core import bioutils
from core import models

def publication_get_form(wiz, form):
    """Constructs the form for publication selection step."""
    user = wiz.request.user
    curator, _ = models.Curator.objects.get_or_create(user=user)
    assigned_pubs = models.Publication.objects.filter(assigned_to=curator,
                                                      curation_complete=False)
    form.fields['pub'].choices = [
        (p.publication_id, '; '.join((p.title, p.authors, p.journal,
                                      'PMID: %s' % p.pmid)))
        for p in assigned_pubs]

    # External curators shouldn't see "This paper contains no data" checkbox
    # on their paper list, assuming that they are submitting on a paper they
    # uploaded and which contains data.
    if curator.curator_type == 'external':
        form.fields['no_data'].initial = False
        form.fields['no_data'].widget = forms.HiddenInput()

    return form

def genome_get_form(wiz, form):
    """Constructs the form for genome and TF selection step."""
    prev_curation = wiz.request.session.get('previous_curation')
    # If selected publication is the one most recently curated, the related
    # curation should be in object prev_curation. Otherwise, it is None.  If so,
    # populate 'Genome' and 'TF information' form fields from the previously
    # submitted curation to make things easier for curator.
    if prev_curation:
        form.initial['TF'] = prev_curation.TF
        sites = prev_curation.site_instances.all()
        if sites:
            genomes = list(set(site.genome.genome_accession for site in sites))
            form.initial['genome_accession'] = sites[0].genome.genome_accession
            for i, genome in enumerate(genomes[1:], start=1):
                form.initial['genome_accession_%d' % i] = genome

        TF_instances = prev_curation.TF_instances.all()
        #form.initial['TF_accession'] = TF_instances[0].protein_accession
        #for i, TF_instance in enumerate(TF_instances):
        #    form.initial['TF_accession_%d' % i] = TF_instance.protein_accession

        form.initial["TF_species"] = prev_curation.TF_species
        form.initial["site_species"] = prev_curation.site_species

        messages.warning(wiz.request, mark_safe("""
        <h4>Warning!</h4> It seems that the paper you selected is
        previously curated. For convenience, fields in this form are
        automatically filled based on the previous curation of the paper. They
        may differ in this curation, so it is best to check that they are
        correct before proceeding to the next step."""))

    # In addition populate two fields on whether the manuscript contains
    # experimental data and promoter information
    pub_id = wiz.request.session['publication']
    p = models.Publication.objects.get(publication_id=pub_id)
    form.fields['contains_promoter_data'].initial = p.contains_promoter_data
    form.fields['contains_expression_data'].initial = p.contains_expression_data

    return form

def techniques_get_form(wiz, form):
    """Construct the form for experiemental techniques step."""
    c = wiz.request.session.get('previously_curated_paper')
    # if selected paper is previously curated, prepopulate experimental
    # techniques
    if c:
        # get all techniques used in this curation
        cur_site_insts = models.Curation_SiteInstance.objects.filter(curation=c)
        techniques = list(set(t.technique_id for csi in cur_site_insts
                              for t in csi.experimental_techniques.all()))
        form.fields['techniques'].initial = techniques

        form.fields['experimental_process'].initial = c.experimental_process
        try:
            external_dbs = models.Curation_ExternalDatabase.objects.filter(curation=c)
            for i,external_db in enumerate(external_dbs):
                form.fields['external_db_type_%d'%i].initial = external_db.external_database.ext_database_id
                form.fields['external_db_accession_%d'%i].initial = external_db.accession_number
        except models.Curation_ExternalDatabase.DoesNotExist:
            pass
        form.fields['forms_complex'].initial = c.forms_complex
        form.fields['complex_notes'].initial = c.complex_notes
    return form

def site_entry_get_form(wiz, form):
    """Construct the form for site entry step."""
    c = wiz.request.session.get('previously_curated_paper')
    # if paper is previously curated, prepopulate fields
    if c:
        # Delete session data, if user change any field and then come back,
        # Store users last entered data, instead of populated data.
        wiz.request.session['previously_curated_paper'] = None

    # populate motifs for the curator to map sites to one of them.
    genomes = session_utils.get(wiz.request.session, 'genomes')
    TF_instances = session_utils.get(wiz.request.session, 'TF_instances')
    csis = models.Curation_SiteInstance.objects.filter(
        site_type='motif_associated',
        site_instance__genome__in=genomes,
        curation__TF_instances=TF_instances)
    # add motifs as site_type choices
    site_type_choices = []
    motif_ids = csis.values_list('motif_id', flat=True).distinct()
    for motif_id in motif_ids:
        sites = [csi.site_instance
                 for csi in csis.filter(motif_id=motif_id)]
        weblogo = mark_safe("<img src='%s'" %
                            bioutils.weblogo_uri(bioutils.run_lasagna(sites)))
        site_type_choices.append((motif_id, weblogo))

    # add non-motif-associated and variable-motif-associated
    site_type_choices.append((max(motif_ids)+1 if motif_ids else 0,
                              'motif-associated (new motif)'))
    site_type_choices.append(('var_motif_associated',
                              'variable motif associated'))
    site_type_choices.append(('non_motif_associated',
                              'non-motif associated'))
    

    form.fields['site_type'].choices = site_type_choices

    # if not high-throughput mode, delete related fields
    if not wiz.request.session.get('high_throughput_curation'):
        del form.fields['peaks']
        del form.fields['assay_conditions']
        del form.fields['method_notes']
        del form.fields['peak_techniques']
    else: # high-throughput mode
        # Populate peak techniques
        techniques = wiz.request.session.get('techniques')
        choices = [(t.technique_id, t.name) for t in techniques]
        form.fields['peak_techniques'].choices = choices


    return form

def site_exact_match_get_form(wiz, form):
    """Show the list of sites and their exact matches in this form.  Each
    reported site and its exact match results are represented as a field. The
    form is generated dynamically as the number of exact matches is not
    fixed."""
    sites = wiz.request.session.get('sites')
    for site in sites:
        label = site.seq
        choices = site.populate_match_choices(True, 'exact_only')
        form.fields[site.key] = forms.ChoiceField(label=label, choices=choices,
                                                  widget=forms.RadioSelect())
        form.fields[site.key].initial = str(choices[0][0])

    return form

def site_soft_match_get_form(wiz, form):
    """For sites that are not matched in the previous step (exact match step),
    show soft-search results."""
    sites = wiz.request.session.get('sites')
    for site in sites:
        # Render all site that are not exactly matched in the previous step
        if not (site.is_matched() and site.get_match().is_exact()):
            label = site.seq
            choices = site.populate_match_choices(True, 'inexact_only')
            form.fields[site.key] = forms.ChoiceField(label=label, choices=choices,
                                                      widget=forms.RadioSelect())
            form.fields[site.key].initial = str(choices[0][0])
    return form

def site_annotation_get_form(wiz, form):
    """Annotation step for site instances."""
    sites = wiz.request.session.get('sites')
    techniques = wiz.request.session.get('techniques')
    for site in sites:
        if site.is_matched() or True:
            i = site.key
            # create dummy field for the label
            label = (site.get_match().pprint() if site.is_matched()
                     else site.pprint())
            form.fields['%d_site' % i] = forms.BooleanField(label=label,
                                                            required=False)
            # create quantitative value field
            if wiz.request.session.get('has_quantitative_data'):
                form.fields['%d_qval' % i] = forms.FloatField(label='Q-val',
                                                              required=False,
                                                              initial=site.qval)
            # create TF type selection field
            form.fields['%d_TF_type' % i] = forms.ChoiceField(label='TF type',
                                                              choices=models.Curation_SiteInstance.TF_TYPE)
            # create TF function selection field
            form.fields['%d_TF_function' % i] = forms.ChoiceField(label='TF function',
                                                                  choices=models.Curation_SiteInstance.TF_FUNCTION)
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
    pid = wiz.request.session.get('publication')
    publication = models.Publication.objects.get(pk=pid)
    # Get all site matches
    sites = wiz.request.session.get('sites')
    for site in sites:
        if site.is_matched():
            i = site.key
            choices = [(g.gene_id, "%s (%s): %s" %
                        (g.locus_tag, g.name, g.description))
                       for g in site.get_match().nearby_genes]
            form.fields[i] = forms.MultipleChoiceField(
                label=site.get_match().pprint(False),
                choices=choices,
                required=False,
                widget=forms.CheckboxSelectMultiple(),
                help_text=site.get_match().match_diagram)

            if not publication.contains_expression_data:
                form.fields[i].widget.attrs['disabled'] = 'disabled'

    return form

def curation_review_get_form(wiz, form):
    """Get curation review form"""
    curator = models.Curator.objects.get(user=wiz.request.user)
    # If curator is not admin, mark curation as requiring validation
    if curator.curator_type == "external":
        form.fields['confidence'].widget = forms.HiddenInput()
    return form
