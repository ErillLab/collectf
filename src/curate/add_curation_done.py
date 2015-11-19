
"""=done= method in Django wizard-view specifies what should happen when the
data for every form is submitted and validated.

Since CollecTF curation form is fairly complex and consists of many steps,
=done= function is split into "sub-done" methods, one for each form-step.
"""

from django.contrib import messages
from django.shortcuts import redirect

from core import models
from . import session_utils
from .forms.add_curation_form import GenomeForm
from .forms.add_curation_form import TechniquesForm
from .forms.add_curation_form import SiteEntryForm
from .forms.add_curation_form import CurationReviewForm


def head(xs):
    """Returns the first element of the list"""
    return xs[0]


def publication_done(wiz):
    """The first step is fairly simple. All the curator does in this step is to
    select the paper to be curated."""
    pid = session_utils.get(wiz.request.session, 'publication')
    publication = models.Publication.objects.get(pk=pid)
    return publication


def genome_done(wiz, form_list):
    """Get genome and TF data from the genome. (2nd step)"""
    form = head([f for f in form_list if type(f) == GenomeForm])

    # Get TF instances
    TF_accession = form.cleaned_data['TF_accession']
    TF_instances = [models.TFInstance.objects.get(
        uniprot_accession=TF_accession)]
    # Extra genome accession numbers (if any)
    for i in xrange(form.NUM_EXTRA_TF_FIELDS):
        tf = form.cleaned_data.get('TF_accession_%d' % i, None)
        if tf:
            TF_instances.append(models.TFInstance.objects.get(
                uniprot_accession=tf))

    assert TF_instances
    return dict(TF=form.cleaned_data['TF'],
                TF_species=form.cleaned_data['TF_species'],
                site_species=form.cleaned_data['site_species'],
                TF_instances=TF_instances)


def techniques_done(wiz, form_list):
    """In this step, the list of used techniques are reported. Since used
    techniques are linked to site instances in the "site-annotation" step,
    there is no need to get any data from this step, except the overall
    description.

    """
    form = head([f for f in form_list if type(f) == TechniquesForm])
    d = dict(experimental_process=form.cleaned_data['experimental_process'],
             forms_complex=form.cleaned_data['forms_complex'],
             complex_notes=form.cleaned_data['complex_notes'])
    for i in xrange(form.NUM_EXTRA_DB_FIELDS):
        d['external_db_type_%d' % i] = form.cleaned_data.get(
            'external_db_type_%d' % i, None)
        d['external_db_accession_%d' % i] = form.cleaned_data.get(
            'external_db_accession_%d' % i, None)
    return d


# The following few steps are related to entering sites and mapping them to
# sequences from the genome. They are carried through the curation using Django
# cache utility. There is no data required from individual steps about site
# instances. All the site-instance data is processed together.
def site_entry_done(wiz, form_list):
    form = head([f for f in form_list if type(f) == SiteEntryForm])
    d = dict(site_type=form.cleaned_data['site_type'])
    # If the curation is high-throughput, capture assay conditions and method
    # notes
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        d['method_notes'] = form.cleaned_data['method_notes']
        d['assay_conditions'] = form.cleaned_data['assay_conditions']
        d['peak_techniques'] = form.cleaned_data['peak_techniques']

    return d


def site_exact_match_done(wiz, form_list):
    pass


def site_soft_match_done(wiz, form_list):
    pass


def site_annotation_done(wiz, form_list):
    pass


def gene_regulation_done(wiz, form_list):
    pass


def curation_review_done(wiz, form_list):
    form = head([f for f in form_list if type(f) == CurationReviewForm])
    return {
        'requires_revision': form.cleaned_data['revision_reasons'],
        'confidence': form.cleaned_data['confidence'],
        'paper_complete': form.cleaned_data['paper_complete'],
        'notes': form.cleaned_data['notes']}


def master_done(wiz, form_list, **kwargs):
    """Done function which calls step-specific done functions and put
    everything together."""
    genome_dict = genome_done(wiz, form_list)  # Genome-step data
    techniques_dict = techniques_done(wiz, form_list)  # Techniques-step data
    site_entry_dict = site_entry_done(wiz, form_list)  # Site-report-step data

    # Curation-review-step data
    curation_review_dict = curation_review_done(wiz, form_list)

    # Create curation object
    curation = create_curation(wiz, genome_dict, techniques_dict,
                               curation_review_dict)
    # Add external db (if any)
    add_external_db(techniques_dict, curation)

    # If it is high-throughput submission add ChIP-info notes
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        add_high_throughput_notes(site_entry_dict, curation)

    # Create matched and non-matched site instances and their related
    # regulation objects
    create_site_instances(wiz, curation, site_entry_dict['site_type'])

    # Mark the paper as complete if so
    paper_complete(wiz, curation_review_dict)

    # Clear session
    session_utils.clear(wiz.request.session)

    # Return success message
    messages.success(wiz.request, "Curation was successfully submitted.")
    return redirect('view_curation', curation.pk)


def is_ncbi_ready(wiz):
    """Check if the curation is ready to submit to NCBI.  If at least 90% of
    reported sites have exact matches in the reference genome, the curation is
    marked as ready for NCBI submission."""
    sites = session_utils.get(wiz.request.session, 'sites')
    exacts = [site
              for site in sites
              if site.is_matched() and site.get_match().is_exact()]
    if len(exacts) < 0.9 * len(sites):
        return False
    return True

def create_curation(wiz, genome_dict, techniques_dict, curation_review_dict):
    """Creates curation object and save it to the database."""
    # Find the curator
    curator = models.Curator.objects.get(user=wiz.request.user)
    # Get publication information
    publication = publication_done(wiz)
    # Create curation
    curation = models.Curation(
        publication=publication,
        TF_species=genome_dict['TF_species'],
        site_species=genome_dict['site_species'],
        TF=genome_dict['TF'],
        experimental_process=techniques_dict['experimental_process'],
        forms_complex=techniques_dict['forms_complex'],
        complex_notes=techniques_dict['complex_notes'],
        requires_revision=curation_review_dict['requires_revision'],
        notes=curation_review_dict['notes'],
        confidence=curation_review_dict['confidence'],
        NCBI_submission_ready=is_ncbi_ready(wiz),
        curator=curator)
    curation.save()

    # If the curation has an associated quantitative data format, add it
    qformat = session_utils.get(wiz.request.session,
                                'quantitative_data_format')
    curation.quantitative_data_format = qformat
    curation.save()

    # Add TF instances
    for TF_instance in genome_dict['TF_instances']:
        curation.TF_instances.add(TF_instance)
    curation.save()

    return curation


def create_site_instances(wiz, curation, site_type):
    """Save Curation_SiteInstance objects."""
    # Get sites from the session data.
    sites = session_utils.get(wiz.request.session, 'sites')
    matched_sites = [site for site in sites if site.is_matched()]
    non_matched_sites = [site for site in sites if not site.is_matched()]
    # Create matched Curation_SiteInstances
    create_matched_site_instances(wiz, curation, matched_sites, site_type)
    # Create non-matched sites
    create_non_matched_site_instances(wiz, curation, non_matched_sites)
    # If there is any peak data, save them as non-motif-associated
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        peaks = session_utils.get(wiz.request.session, 'peaks')
        create_peak_site_instances(wiz, curation, peaks)


def get_site_type_and_motif_id(site_type):
    """The site_type field on curation wizard is used to figure out site type
    and motif id.

    If the site_type is selected as var_motif_associated or
    non_motif_associated, the sites are marked as variable motif associated and
    non-motif associated, respectively. If a motif is selected, sites are
    marked as motif-associated and the motif id is set accordingly.
    """
    if site_type in ['var_motif_associated', 'non_motif_associated']:
        return site_type, -1    # motif id is not important any more

    return 'motif_associated', int(site_type)


def create_matched_site_instances(wiz, curation, sites, site_type):
    """Create SiteInstance objects and CurationSiteInstance objects to link
    them to the Curation.
    """
    for site in sites:
        # Create SiteInstance objects (or get if they are already created).
        assert site.is_matched(), "Site is not matched."
        match = site.get_match()
        s, _ = models.SiteInstance.objects.get_or_create(
            genome=match.genome,
            start=match.start,
            end=match.end,
            strand=match.strand,
            _seq=match.seq)

        # Create Curation_SiteInstance object
        stype, motif_id = get_site_type_and_motif_id(site_type)
        cs = models.Curation_SiteInstance(
            curation=curation,
            site_instance=s,
            annotated_seq=match.reported_seq,
            site_type=stype,
            motif_id=motif_id,
            TF_function=site.TF_function,
            TF_type=site.TF_type,
            quantitative_value=site.qval)
        cs.save()

        # For each site instance, add experimental technique information
        for t in site.techniques:
            cs.experimental_techniques.add(t)
        cs.save()

        # For site instance, add regulation information
        for gene in match.nearby_genes:
            # If a nearby gene is marked as having experimental evidence, it is
            # saved as "experimentally verified", otherwise a nearby gene is
            # saved as "inferred" regulation.
            evidence = ('exp_verified' if gene in match.regulated_genes
                        else 'inferred')
            models.Regulation(
                curation_site_instance=cs,
                gene=gene,
                evidence_type=evidence).save()


def create_non_matched_site_instances(wiz, curation, sites):
    """Create not matched site instances as NotAnnotatedSiteInstance
    objects."""
    for site in sites:
        assert not site.is_matched(), "Site is matched."
        nms = models.NotAnnotatedSiteInstance(
            sequence=site.seq,
            curation=curation,
            TF_function=site.TF_function,
            TF_type=site.TF_type,
            quantitative_value=site.qval)
        nms.save()
        # For each site instance, add experimental techniques
        for tech in site.techniques:
            nms.experimental_techniques.add(tech)
        nms.save()


def create_peak_site_instances(wiz, curation, peaks):
    """Create peaks as non-motif-associated site instances"""
    matched_peaks = [p for p in peaks if p.is_matched()]
    non_matched_peaks = [p for p in peaks if not p.is_matched()]
    create_matched_peaks(wiz, curation, matched_peaks)
    create_non_matched_peaks(wiz, curation, non_matched_peaks)


def create_matched_peaks(wiz, curation, peaks):
    for peak in peaks:
        assert peak.is_matched()
        match = peak.get_match()
        s, _ = models.SiteInstance.objects.get_or_create(
            genome=match.genome,
            start=match.start,
            end=match.end,
            strand=match.strand,
            _seq=match.seq)

        # Create the curation_site_instance object
        cs = models.Curation_SiteInstance(
            curation=curation,
            site_instance=s,
            annotated_seq=match.reported_seq,
            site_type='non_motif_associated',
            motif_id=-1,
            TF_function='N/A',
            TF_type='N/A',
            quantitative_value=peak.qval,
            is_high_throughput=True)
        cs.save()

        # For each peak instance, add the experimental technique information
        for t in peak.techniques:
            cs.experimental_techniques.add(t)
        cs.save()

        # For site instance, add regulation information
        for gene in match.nearby_genes:
            # If a nearby gene is marked as having experimental evidence, it is
            # saved as "experimentally verified", otherwise a nearby gene is
            # saved as "inferred" regulation.
            evidence = ('exp_verified' if gene in match.regulated_genes
                        else 'inferred')
            models.Regulation(curation_site_instance=cs,
                              gene=gene,
                              evidence_type=evidence).save()


def create_non_matched_peaks(wiz, curation, peaks):
    for peak in peaks:
        assert not peak.is_matched()
        nms = models.NotAnnotatedSiteInstance(
            curation=curation,
            sequence=peak.seq,
            TF_function='N/A',
            TF_type='N/A',
            quantitative_value=peak.qval,
            is_high_throughput=True)
        nms.save()
        for t in peak.techniques:
            nms.experimental_techniques.add(t)
        nms.save()


def add_external_db(techniques_dict, curation):
    """Add external DB references, if any."""
    for i in range(5):
        external_db_type = techniques_dict.get('external_db_type_%d' % i, None)
        if external_db_type and external_db_type != 'None':
            external_db_type = models.ExternalDatabase.objects.get(
                ext_database_id=techniques_dict['external_db_type_%d' % i])
            curation_ext_ref = models.Curation_ExternalDatabase(
                curation=curation,
                external_database=external_db_type,
                accession_number=techniques_dict[
                    'external_db_accession_%d' % i])
            curation_ext_ref.save()


def add_high_throughput_notes(site_entry_dict, curation):
    """If the curation is high-throughput submission, save assay conditions and
    method notes in the database and link it to the curation
    """
    assay_conditions = site_entry_dict['assay_conditions']
    method_notes = site_entry_dict['method_notes']
    chip_info = models.ChipInfo(assay_conditions=assay_conditions,
                                method_notes=method_notes)
    chip_info.save()
    curation.chip_info = chip_info
    curation.save()


def paper_complete(wiz, curation_review_dict):
    """Mark the paper complete if so."""
    # Get publication information
    publication = publication_done(wiz)
    if curation_review_dict['paper_complete']:
        publication.curation_complete = True
        publication.save()
