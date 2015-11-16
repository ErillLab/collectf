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
