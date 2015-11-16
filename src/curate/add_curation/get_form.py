from django import forms
from django.contrib import messages
from django.utils.safestring import mark_safe

from core import models


def publication_get_form(wiz, form):
    """Constructs the form for publication selection step."""
    user = wiz.request.user
    curator, _ = models.Curator.objects.get_or_create(user=user)
    assigned_pubs = models.Publication.objects.filter(assigned_to=curator,
                                                      curation_complete=False)
    form.fields['publication'].choices = [
        (pub.publication_id, '; '.join([pub.title, pub.authors, pub.journal]))
        for pub in assigned_pubs]

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
            form.initial['genome_accession'] = sites[0].genome.genome_accession

        TF_instances = prev_curation.TF_instances.all()
        form.initial['TF_accession'] = TF_instances[0].protein_accession
        form.initial['TF_species'] = prev_curation.TF_species
        form.initial['site_species'] = prev_curation.site_species

        messages.warning(wiz.request, mark_safe("""
        <h4>Warning!</h4> It seems that the paper you selected is
        previously curated. For convenience, fields in this form are
        automatically filled based on the previous curation of the paper. They
        may differ in this curation, so it is best to check that they are
        correct before proceeding to the next step."""))

    # Populate two additional fields on whether the manuscript contains
    # experimental data and promoter information
    pid = wiz.request.session.get('publication')
    pub = models.Publication.objects.get(publication_id=pid)
    form.fields['contains_promoter_data'].initial = pub.contains_promoter_data
    form.fields['contains_expression_data'].initial = pub.contains_expression_data

    return form
