from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from collectfapp.curationview import CurationWizard
from collectfapp.curationform import *
import models
import sutils
import views

def init_publication(curation):
    """Return publication data for new curation obj."""
    return dict(pub=curation.publication.publication_id)

def get_genome(curation):
    # get genome
    site_instances = curation.site_instances.all()
    # make sure there is >1 site instance for curation. No other way to access genome
    # Wed Feb 27 18:07:02 2013 => Actually, it is possible for a curation to not to
    # have any site. In this case, there is no way to access the genome accession
    # number reported in the original curation.
    if site_instances:
        genome_accession = site_instances[0].genome.genome_accession
    else:
        genome_accession = ""

    return genome_accession

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
    return [str(t.technique_id) for t in ts]

def init_techniques_form(curation):
    """Return techniques form data from old curation for new one."""
    try:
        external_db = models.Curation_ExternalDatabase.objects.get(curation=curation)
    except Curation_ExternalDatabase.DoesNotExist:
        external_db = None
        
    return dict(techniques = get_used_techniques(curation),
                experimental_process = curation.experimental_process,
                forms_complex = curation.forms_complex,
                complex_notes = curation.complex_notes,
                external_db_type = external_db.external_database.ext_database_id if external_db else None,
                external_db_accession = external_db.accession_number if external_db else "")
        

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
    messages.warning(request, 'Edit curation is temporarily unavailable.')
    return HttpResponseRedirect(reverse(views.home))
    # get curation
    old_curation = models.Curation.objects.get(curation_id=cid)
    # get initial data for new curation form
    initial = {'0': init_publication(old_curation),
               '1': init_genome_form(old_curation),
               '2': init_techniques_form(old_curation),
               '3': init_site_report_form(old_curation),
               '9': init_curation_review_form(old_curation),
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

    # When doing a new curation, the first form is publication step. To help the
    # curator, if the publication is previously curated, some fields are
    # pre-populated, such as TF-type, TF-function, TF-species, etc. (assuming they
    # are same for all curations that belong to one publication). We shouldn't use
    # this feature, for edit_curation feature as the fields are populated with
    # curation data being edited.
    sutils.sput(request.session, "previously_curated_paper", None)

    wiz = CurationWizard.as_view([PublicationForm,
                                  GenomeForm,
                                  TechniquesForm,
                                  SiteReportForm,
                                  SiteExactMatchForm,
                                  SiteSoftMatchForm,
                                  SiteQuantitativeDataForm,
                                  SiteRegulationForm,
                                  CurationReviewForm],
                                 initial_dict=initial,
                                 condition_dict={'0': False,
                                                 '4': exact_site_match_form_condition,
                                                 '5': inexact_site_match_form_condition,
                                                 '6': site_quantitative_data_form_condition})
                                 
    return wiz(request)
    
def exact_site_match_form_condition(wizard):
    return True

def inexact_site_match_form_condition(wizard):
    return sutils.sget(wizard.request.session, 'soft_site_match_choices')

def site_quantitative_data_form_condition(wizard):
    return sutils.sget(wizard.request.session, 'has_quantitative_data')
