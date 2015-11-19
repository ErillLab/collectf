"""This file contains process functions for curation-wizard steps. After
submission of each step, the form is validated. The clean and validate data is
passed to process function of the associated step for post-processing.
"""

from core import models

from . import site_entry
from . import session_utils


def publication_process(wiz, form):
    """Post-processes paper selection step."""
    pubid = form.cleaned_data['pub']
    session_utils.put(wiz.request.session, 'publication', int(pubid))
    if form.cleaned_data['no_data']:
        # mark paper as having no data
        paper = models.Publication.objects.get(publication_id=pubid)
        note = " \nPaper has no TF-binding site data."
        paper.submission_notes += note
        paper.curation_complete = True
        paper.save()
        session_utils.put(wiz.request.session, 'paper_contains_no_data', True)
        return

    # If paper is previously curated, populate genome and TF information form
    # search DB if there is any curation object belonging to this publication
    p = models.Publication.objects.get(publication_id=pubid)
    # Check if the publication is previously curated
    cs = models.Curation.objects.filter(publication=p)
    if cs.count() >= 1:
        session_utils.put(wiz.request.session, 'previous_curation', cs[0])
    else:
        session_utils.remove(wiz.request.session, 'previous_curation')


def genome_process(wiz, form):
    """Post-processes genome and TF selection step."""
    # In form validation step, genome is searched in the database, and if not
    # found, it is inserted into db. So, at this point, it is guaranteed that
    # genome with id <genome_accession> should be in DB.
    genome_accession = form.cleaned_data['genome_accession']
    genomes = [models.Genome.objects.get(genome_accession=genome_accession)]
    # Extra genome accession numbers (if any)
    for i in xrange(form.NUM_EXTRA_GENOME_FIELDS):
        g = form.cleaned_data.get('genome_accession_%d' % i, None)
        if g:
            genomes.append(models.Genome.objects.get(genome_accession=g))

    session_utils.put(wiz.request.session, 'genomes', genomes)

    protein_accessions = [form.cleaned_data['TF_accession']]
    # Extra TF accessions (if any)
    for i in xrange(form.NUM_EXTRA_TF_FIELDS):
        t = form.cleaned_data.get('TF_accession_%d' % i, None)
        if t:
            protein_accessions.append(t)
    TF_instances = models.TFInstance.objects.filter(
        uniprot_accession__in=protein_accessions)

    session_utils.put(wiz.request.session, 'TF_instances', TF_instances)
    session_utils.put(wiz.request.session, 'site_species',
                      form.cleaned_data['site_species'])
    session_utils.put(wiz.request.session, 'TF_species',
                      form.cleaned_data['TF_species'])

    # Set manuscript-related fields (contains_experimental_data and
    # contains_promoter_data). Actually, these fields are defined during adding
    # publication process, but the user is given a chance to edit these fields
    # during the curation process.
    pubid = session_utils.get(wiz.request.session, 'publication')
    p = models.Publication.objects.get(publication_id=pubid)
    p.contains_promoter_data = form.cleaned_data["contains_promoter_data"]
    p.contains_expression_data = form.cleaned_data["contains_expression_data"]
    p.save()


def techniques_process(wiz, form):
    """Post-process experimental techniques step."""
    techniques = models.ExperimentalTechnique.objects.filter(
        pk__in=form.cleaned_data['techniques'])
    # save selected techniques (to be used in site-annotation step)
    session_utils.put(wiz.request.session, 'techniques', techniques)


def site_entry_process(wiz, form):
    """Post process site entry step"""
    genomes = session_utils.get(wiz.request.session, 'genomes')
    sites = site_entry.parse_input(form.cleaned_data['sites'].strip())

    # Find exact matches
    for site in sites:
        site.search_exact_match(genomes)

    # If any site has quantitative data, mark the curation
    has_qdata = any(site.qval for site in sites)

    # If high-throughput get peak data to save them as non-motif-associated
    # data
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        peaks = site_entry.parse_input(form.cleaned_data['peaks'].strip())
        techniques = models.ExperimentalTechnique.objects.filter(
            pk__in=form.cleaned_data['peak_techniques'])
        for peak in peaks:
            peak.search_exact_match(genomes)
            # if there is any match, select the first one by default
            if peak.get_exact_matches():
                peak.set_exact_match(0)
            # for each peak add the technique
            peak.clear_techniques()
            for t in techniques:
                peak.add_technique(t)

        if any(peak.qval for peak in peaks):
            has_qdata = True

        session_utils.put(wiz.request.session, 'peaks', peaks)

    # Save the list of sites
    session_utils.put(wiz.request.session, 'sites', sites)
    # Save the type of lists
    session_utils.put(wiz.request.session, 'site_type',
                      form.cleaned_data['site_type'])
    # Save whether curation has quantitative data
    session_utils.put(wiz.request.session, 'has_quantitative_data', has_qdata)
    # If any quantitative data format save it
    qval = form.cleaned_data.get('quantitative_data_format', None)
    session_utils.put(wiz.request.session, 'quantitative_data_format', qval)


def site_exact_match_process(wiz, form):
    """Post process for site exact match step. Identify the sites that are
    matched by one of their possible matches and for the rest, perform a soft
    search which allows some substitutions when searching the sequence in the
    genome"""
    genomes = session_utils.get(wiz.request.session, 'genomes')
    sites = session_utils.get(wiz.request.session, 'sites')
    for site_id, match_id in form.cleaned_data.items():
        site = [site for site in sites if site.key == site_id][0]
        if match_id != 'None':  # means this site is matched
            site.set_exact_match(match_id)
        else:                   # not matched, perform soft search
            site.search_soft_match(genomes)

    # If high-throughput submission, try to match quantitative values in peak
    # data to sites
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        peaks = session_utils.get(wiz.request.session, 'peaks')
        for site in sites:
            site.match_peak_data(peaks)

    # save the list of sites
    session_utils.put(wiz.request.session, 'sites', sites)


def site_soft_match_process(wiz, form):
    """In this form, soft-search results are processed, user matched all of
    sites to any appropriate sequence found in the genome. Some sites might be
    unmatched if there is no any good candidate in search results."""
    sites = session_utils.get(wiz.request.session, 'sites')
    for site_id, match_id in form.cleaned_data.items():
        site = [site for site in sites if site.key == site_id][0]
        if match_id != 'None':  # means this site is matched
            site.set_soft_match(match_id)

    # If high-throughput submission, try to match quantitative values in peak
    # data to sites
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        peaks = session_utils.get(wiz.request.session, 'peaks')
        for site in sites:
            site.match_peak_data(peaks)

    # save the list of sites
    session_utils.put(wiz.request.session, 'sites', sites)


def site_annotation_process(wiz, form):
    """In this form, each individual sites (only matched ones) are further
    annotated. For each site, TF-function (activator, repressor, etc.),
    quantitative value (if any) and experimental techniques used to identify
    site are requested from the curator."""
    sites = session_utils.get(wiz.request.session, 'sites')
    techniques = session_utils.get(wiz.request.session, 'techniques')
    has_qdata = session_utils.get(wiz.request.session, 'has_quantitative_data')
    cd = form.cleaned_data
    for site in sites:
        i = site.key

        # Make sure all matched sites are  in the annotation form
        assert '%d_site' % i in cd, "Inconsistent site annotation form"

        # Quantitative value
        if has_qdata:
            q = cd['%d_qval' % i]
            site.set_qval(float(q) if q else None)
        # TF function
        site.set_TF_function(cd['%d_TF_function' % i])
        # TF_type
        site.set_TF_type(cd['%d_TF_type' % i])
        # Experimental techniques
        site.clear_techniques()
        for j, t in enumerate(techniques):
            # Add technique to the site if it is checked.
            if cd['%d_technique_%d' % (i, j)]:
                site.add_technique(t)
    # Save sites again
    session_utils.put(wiz.request.session, 'sites', sites)


def gene_regulation_process(wiz, form):
    """Prior to this form all matches sites are processed to identify nearby
    genes. In this step, nearby genes are displayed and asked to the curator to
    select genes that have experimental evidence of regulation in the curated
    paper."""
    cd = form.cleaned_data
    sites = session_utils.get(wiz.request.session, 'sites')
    for site in sites:
        i = site.key
        if not site.is_matched():
            continue

        # Make sure this site is in annotation form
        assert i in cd, "Inconsistent gene regulation form."
        match = site.get_match()
        match.set_regulated_genes(models.Gene.objects.filter(pk__in=cd[i]))


def review_curation_process(wiz, form):
    pass
