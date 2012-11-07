from django.contrib.auth.decorators import login_required
from collectfapp.curationview import CurationWizard
from collectfapp.curationform import *
import models
import sutils

def init_publication(curation):
    """Return publication data for new curation obj."""
    return dict(pub=curation.publication.publication_id)

def get_genome(curation):
    # get genome
    site_instances = curation.site_instances.all()
    # make sure there is >1 site instance for curation. No other way to access genome
    assert site_instances
    return site_instances[0].genome.genome_accession

def init_genome_form(curation):
    """Return genome form data from old curation for new one"""
    return dict(genome_accession = get_genome(curation),
                TF_accession = curation.TF_instance.protein_accession,
                TF = curation.TF.TF_id if curation.TF else None,
                TF_type = curation.TF_type,
                TF_function = curation.TF_function,
                TF_species = curation.TF_species,
                site_species = curation.site_species)

def get_used_techniques(curation):
    # return curation techniques
    ts = curation.experimental_techniques.all()
    return [t.technique_id for t in ts]

def init_techniques_form(curation):
    """Return techniques form data from old curation for new one."""
    return dict(techniques = get_used_techniques(curation),
                experimental_process = curation.experimental_process,
                forms_complex = curation.forms_complex,
                complex_notes = curation.complex_notes)

def init_site_report_form(curation):
    """Return reported sites from old curation"""
    # get matched sites
    seqs = [s.annotated_seq for s in curation.curation_siteinstance_set.all()]
    # and not matched/annotated sites
    not_sites = models.NotAnnotatedSiteInstance.objects.filter(curation=curation)
    not_seqs = [s.sequence for s in not_sites]
    return dict(sites = '\n'.join(seqs + not_seqs))

def init_curation_review_form(curation):
    """Return curation review data from the old curation"""
    return dict(revision_reasons = curation.requires_revision,
                confidence = curation.confidence,
                paper_complete = curation.publication.curation_complete,
                NCBI_submission_ready = curation.NCBI_submission_ready,
                notes = curation.notes)


@login_required
def edit_curation(request, cid):
    """Handler function for curation editing.
    - Get curation being edited.
    - Create new form wizard with all initial data from old curation (except site
    matches data)
    - When new form is submitted, insert new curation instance to DB and remove the
    old one.
    """
    # get curation
    old_curation = models.Curation.objects.get(curation_id=cid)
    # get initial data for new curation form
    initial = {'0': init_publication(old_curation),
               '1': init_genome_form(old_curation),
               '2': init_techniques_form(old_curation),
               '3': init_site_report_form(old_curation),
               '7': init_curation_review_form(old_curation),
              }
    
    # since if paper complete, it will not be displayed in the first form as
    # available. Since we are using CurationForm for edit curation too, change it to
    # false.
    #old_curation.publication.curation_complete = False
    #old_curation.publication.save()

    # tell form wizard that this is an edit to an existing curation
    sutils.sput(request.session, 'old_curation', old_curation)

    # save which publication we are about to edit
    sutils.sput(request.session, 'publication', old_curation.publication.publication_id)

    wiz = CurationWizard.as_view([PublicationForm,
                                  GenomeForm,
                                  TechniquesForm,
                                  SiteReportForm,
                                  SiteExactMatchForm,
                                  SiteSoftMatchForm,
                                  SiteRegulationForm,
                                  CurationReviewForm],
                                 initial_dict=initial,
                                 condition_dict={'0': False}) # Don't show publication form
                                 
    return wiz(request)
    

    
    
