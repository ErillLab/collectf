from core import models


def publication_process_step(wiz, form):
    """Post-processes paper selection step."""

    publication_id = int(form.cleaned_data['publication'])
    wiz.request.session['publication'] = publication_id
    paper = models.Publication.objects.get(publication_id=publication_id)

    if form.cleaned_data.get('no_data'):
        # Mark paper as having no data
        note = " \nPaper has no TF-binding site data."
        paper.submission_notes += note
        paper.curation_complete = True
        paper.save()
        wiz.request.session['paper_contains_no_data'] = True
        return

    # If paper is previously curated, populate genome and TF information form
    # search DB if there is any curation object belonging to this publication

    # Check if the publication is previously curated.
    curations = models.Curation.objects.filter(publication=paper)
    if curations:
        wiz.request.session['previous_curation'] = curations.all()[0]
    else:
        wiz.request.session['previous_curation'] = None


def genome_process_step(wiz, form):
    """Post-process genome and TF selection step."""
    
    # Genome is searched in database during form validation, and if not found,
    # it is saved into the database.. So, at this point, it is guaranteed that
    # genome with id <genome_accession> should be in db.
    genome_accessions = [form.cleaned_data['genome_accession']]
    # Extra genome accession numbers (if any)
    for i in range(form.NUM_EXTRA_GENOME_FIELDS):
        acc = form.cleaned_data.get('genome_accession_%d' % i, None)
        if acc:
            genome_accessions.append(acc)
    
    # Store genome in session data
    genomes = models.Genome.objects.filter(
        genome_accession__in=genome_accessions)
    wiz.request.session['genomes'] = genomes

    TF_accessions = [form.cleaned_data['TF_accession']]
    # Extra TF accessions (if any)
    for i in range(form.NUM_EXTRA_TF_FIELDS):
        acc = form.cleaned_data.get('TF_accession_%d' % i, None)
        if acc:
            TF_accessions.append(t)
            
    # Store TF accessions in session data
    TF_instances = models.TFInstance.objects.filter(
        uniprot_accession__in=TF_accessions)
    wiz.request.session['TF_instances'] = TF_instances

    wiz.request.session['site_species'] = form.cleaned_data['site_species']
    wiz.request.session['TF_species'] = form.cleaned_data['TF_species']

    # Set manuscript-related fields (contains_experimental_data and
    # contains_promoter_data). These fields are also shown in add-publication
    # form, but the user is given a chance to edit these fields during
    # curation.
    pubid = wiz.request.session['publication']
    pub = models.Publication.objects.get(publication_id=pubid)
    pub.contains_promoter_data = form.cleaned_data['contains_promoter_data']
    pub.contains_expression_data = form.cleaned_data['contains_expression_data']
    pub.save()


def techniques_process_step(wiz, form):
    """Post-processes experimental techniques step."""
    
    techniques = models.ExperimentalTechnique.objects.filter(
        pk__in=form.cleaned_data['techniques'])
    # Save selected techniques (to be used in site-annotation step)
    wiz.request.session['techniques'] = techniques
