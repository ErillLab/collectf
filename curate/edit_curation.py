"""View for editing curations."""

from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.contrib import messages
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from add_curation import CurationWizard
from forms.add_curation_form import *
from base import models
import session_utils

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
    genome_accessions = []

    if site_instances:
        genome_accessions = list(set(si.genome.genome_accession for si in site_instances))
    else:
        genome_accessions = []

    return genome_accession


def init_genome_form(curation):
    """Return genome form data from old curation for new one"""
    d = dict(TF = curation.TF.TF_id if curation.TF else None,
             TF_type = curation.TF_type,
             TF_function = curation.TF_function,
             TF_species = curation.TF_species,
             site_species = curation.site_species)
    # Add genome accession numbers
    genomes = get_genome(curation)
    if len(genomes) > 0:
        d['genome_accession'] = genomes[0]
    if len(genomes) == 2:
        d['genome_accession_1'] = genomes[1]
    if len(genomes) == 3:
        d['genome_accession_2'] = genomes[2]
    # Add TF accession numbers
    d['TF_accession'] = curation.site_instances.all()[0].protein_accession
    if len(curation.site_instances.all()) == 2:
        d['TF_accession_1'] = curation.site_instances.all()[1].protein_accession
    if len(curation.site_instances.all()) == 3:
        d['TF_accession_2'] = curation.site_instances.all()[2].protein_accession
    return d

def get_used_techniques(curation):
    """Get all techniques that are used to identify binding sites."""
    # return curation techniques
    cur_site_insts = models.Curation_SiteInstance.objects.filter(curation=curation)
    ts = [t for t in csi.experimental_techniques.all() for csi in cur_site_insts]
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

@user_passes_test(lambda u: u.is_superuser)
def edit_curation(request, cid):
    """Handler function for curation editing.
    - Get curation being edited.
    - Create new form wizard with all initial data from old curation (except site
    matches data)
    - When new form is submitted, insert new curation instance to DB and remove the
    old one.
    """
    # Get curation
    old_curation = models.Curation.objects.get(curation_id=cid)
    # Get initial data for new curation form
    initial = {'0': init_publication(old_curation),
               '1': init_genome_form(old_curation),
               '2': init_techniques_form(old_curation),
               '3': init_site_report_form(old_curation),
               '8': init_curation_review_form(old_curation),
              }
    
    # Tell form wizard that this is an edit to an existing curation
    session_utils.put(request.session, 'old_curation', old_curation)

    # Save which publication we are about to edit
    session_utils.put(request.session, 'publication', old_curation.publication.publication_id)

    # When doing a new curation, the first form is publication step. To help the
    # curator, if the publication is previously curated, some fields are
    # pre-populated, such as TF-type, TF-function, TF-species, etc. (assuming they
    # are same for all curations that belong to one publication). We shouldn't use
    # this feature, for edit_curation feature as the fields are populated with
    # curation data being edited.
    sesion_utils.put(request.session, "previously_curated_paper", None)

    wiz = CurationWizard.as_view([PublicationForm,
                                  GenomeForm,
                                  TechniquesForm,
                                  SiteEntryForm,
                                  SiteExactMatchForm,
                                  SiteSoftMatchForm,
                                  SiteAnnotationForm,
                                  GeneRegulationForm,
                                  CurationReviewForm],
                                 initial_dict=initial,
                                 condition_dict={'0': False})
                                 
    return wiz(request)
