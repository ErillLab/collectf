"""=done= method in Django wizard-view specifies what should happen when the
data for every form is submitted and validated.

Since CollecTF curation form is fairly complex and consists of many steps,
=done= function is split into "sub-done" methods, one for each form-step.
"""
import session_utils
from forms.add_curation_form import *
import base
from base import models
from django.contrib import messages
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse

import browse.view_curation

head = lambda x: x[0]

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
    TF_instances = [models.TFInstance.objects.get(protein_accession=TF_accession)]
    # Extra genome accession numbers (if any)
    for i in xrange(1, settings.NUMBER_OF_TF_ACCESSION_FIELDS):
        tf = form.cleaned_data.get('TF_accession_%d' % i, None)
        if tf: TF_instances.append(models.TFInstance.objects.get(protein_accession=tf))
            
    assert TF_instances
    return dict(
        TF = form.cleaned_data['TF'],
        TF_species = form.cleaned_data['TF_species'],
        site_species = form.cleaned_data['site_species'],
        TF_instances = TF_instances,
    )

def techniques_done(wiz, form_list):
    """In this step, the list of used techniques are reported. Since used
    techniques are linked to site instances in the "site-annotation" step, there
    is no need to get any data from this step, except the overall description."""
    form = head([f for f in form_list if type(f) == TechniquesForm])
    d = dict(
        experimental_process = form.cleaned_data['experimental_process'],
        forms_complex = form.cleaned_data['forms_complex'],
        complex_notes = form.cleaned_data['complex_notes'],
    )
    for i in xrange(settings.NUMBER_OF_EXTERNAL_DATABASE_FIELDS):
        d['external_db_type_%d' % i] = form.cleaned_data.get('external_db_type_%d' % i, None)
        d['external_db_accession_%d' % i] = form.cleaned_data.get('external_db_accession_%d' % i, None)
        
    return d

"""The following few steps are related to entering sites and mapping them to
sequences from the genome. They are carried through the curation using Django
cache utility. There is no data required from individual steps about site
instances. All the site-instance data is processed together."""

def site_entry_done(wiz, form_list):
    form = head([f for f in form_list if type(f) == SiteEntryForm])
    d = dict(
        site_type = form.cleaned_data['site_type'],
    )
    # If the curation is high-throughput, capture assay conditions and method notes
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        d['method_notes'] = form.cleaned_data['method_notes']
        d['assay_conditions'] = form.cleaned_data['assay_conditions']

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
    return dict(
        requires_revision = form.cleaned_data['revision_reasons'],
        confidence = form.cleaned_data['confidence'],
        NCBI_submission_ready = form.cleaned_data['NCBI_submission_ready'],
        paper_complete = form.cleaned_data['paper_complete'],
        notes = form.cleaned_data['notes'],
    )

def master_done(wiz, form_list, **kwargs):
    """Done function which calls step-specific done functions and put everything
    together."""
    # Check if it is an edit_curation wizard
    form_list = edit_curation_check(wiz, form_list)
    # Get data across forms
    
    genome_dict = genome_done(wiz, form_list) # Genome-step data
    techniques_dict = techniques_done(wiz, form_list) # Techniques-step data
    site_entry_dict = site_entry_done(wiz, form_list)    # Site-report-step data
    curation_review_dict = curation_review_done(wiz, form_list) # Curation-review-step data
    
    # Create curation object
    curation = create_curation(wiz, genome_dict, techniques_dict, curation_review_dict)
    # Add external db (if any)
    add_external_db(techniques_dict, curation)
    # Create matched and non-matched site instances and their related regulation objects
    create_site_instances(wiz, curation, site_entry_dict['site_type'])
    # Mark the paper as complete if so
    paper_complete(wiz, curation_review_dict)
    # clear session
    session_utils.clear(wiz.request.session)
    # Return success message
    messages.success(wiz.request, "Curation was successfully submitted.")
    return HttpResponseRedirect(reverse(browse.view_curation.view_curation, kwargs={'cid': curation.pk}))

def create_curation(wiz, genome_dict, techniques_dict, curation_review_dict):
    """Create curation object and save it to the database."""
    # Find the curator
    curator = models.Curator.objects.get(user=wiz.request.user)
    # Get publication information
    publication = publication_done(wiz)
    # Create curation
    curation = models.Curation(publication=          publication,
                               TF_species=           genome_dict['TF_species'],
                               site_species=         genome_dict['site_species'],
                               TF=                   genome_dict['TF'],
                               experimental_process= techniques_dict['experimental_process'],
                               forms_complex=        techniques_dict['forms_complex'],
                               complex_notes=        techniques_dict['complex_notes'],
                               requires_revision=    curation_review_dict['requires_revision'],
                               notes=                curation_review_dict['notes'],
                               confidence=           curation_review_dict['confidence'],
                               NCBI_submission_ready=curation_review_dict['NCBI_submission_ready'],
                               curator=              curator)
    curation.save()

    # If the curation has an associated quantitative data format, add it
    qformat = session_utils.get(wiz.request.session, 'quantitative_data_format')
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
    create_non_matched_site_instances(wiz, curation, non_matched_sites, site_type)
    # If there is any peak data, save them as non-motif-associated
    if session_utils.get(wiz.request.session, 'high_throughput_curation'):
        peaks = session_utils.get(wiz.request.session, 'peaks')
        create_peak_site_instances(wiz, curation, peaks)

def create_matched_site_instances(wiz, curation, sites, site_type):
    """Create SiteInstance objects and CurationSiteInstance objects to link them
    to the Curation."""
    for site in sites:
        # Create SiteInstance objects (or get if they are already created).
        assert site.is_matched(), "Site is not matched."
        match = site.get_match()
        s,_ = models.SiteInstance.objects.get_or_create(genome=match.genome,
                                                        start=match.start,
                                                        end=match.end,
                                                        strand=match.strand,
                                                        _seq=match.seq)
        # Create Curation_SiteInstance object
        cs = models.Curation_SiteInstance(curation=curation, 
                                          site_instance=s,
                                          annotated_seq=match.reported_seq,
                                          site_type=site_type,
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
            evidence = 'exp_verified' if gene in match.regulated_genes else 'inferred'
            models.Regulation(curation_site_instance=cs, gene=gene, evidence_type=evidence).save()

def create_non_matched_site_instances(wiz, curation, sites, site_type):
    """Create not matched site instances as NotAnnotatedSiteInstance objects."""
    for site in sites:
        assert not site.is_matched(), "Site is matched."
        models.NotAnnotatedSiteInstance(sequence=site.seq, curation=curation).save()

def create_peak_site_instances(wiz, curation, peaks):
    """Create peaks as non-motif-associated site instances"""
    for peak in peaks:
        if not peak.is_matched(): continue
        match = peak.get_match()
        s,_ = models.SiteInstance.objects.get_or_create(genome=match.genome,
                                                        start=match.start,
                                                        end=match.end,
                                                        strand=match.strand,
                                                        _seq=match.seq)
        # Create the curation_site_instance object
        cs = models.Curation_SiteInstance(curation=curation,
                                          site_instance=s,
                                          annotated_seq=match.reported_seq,
                                          site_type='non_motif_associated',
                                          TF_function='N/A',
                                          TF_type='N/A',
                                          quantitative_value=peak.qval)
        cs.save()

def edit_curation_check(wiz, form_list):
    """If curation view is used for revising a curation, the first step
    (publication selection) is hidden in the form-wizard and in the
    form_list. However, in that case, the first element (i.e. form_list[0]) is
    not the publication step anymore. The quick solution is to insert a fake
    element to the form_list."""
    if session_utils.has(wiz.request.session, 'old_curation'):
        # This is a revision for an existing curation
        form_list.insert(0, None)
    return form_list

def add_external_db(techniques_dict, curation):
    """Add external DB references, if any."""
    for i in xrange(settings.NUMBER_OF_EXTERNAL_DATABASE_FIELDS):
        external_db_type = techniques_dict.get('external_db_type_%d' % i, None)
        if external_db_type and external_db_type != 'None':
            external_db_type = models.ExternalDatabase.objects.get(
                ext_database_id=techniques_dict['external_db_type_%d' % i])
            curation_ext_ref = models.Curation_ExternalDatabase(
                curation=curation, external_database=external_db_type,
                accession_number=techniques_dict['external_db_accession_%d' % i])
            curation_ext_ref.save()

def add_high_throughput_notes(site_entry_dict, curation):
    """If the curation is high-throughput submission, save assay conditions and
    method notes in the database and link it to the curation"""
    assay_conditions = site_entry_dict['assay_conditions']
    method_notes = site_entry_dict['method_notes']
    chip_info = models.ChipInfo(assay_conditions=assay_conditions, method_notes=method_notes)
    chip_info.save()
    curation.chip_info = chip_info
    curation.save()
def paper_complete(wiz, curation_review_dict):
    """Mark the paper complete if so"""
    # Get publication information
    publication = publication_done(wiz)
    if curation_review_dict['paper_complete']:
        publication.curation_complete = True
        publication.save()
