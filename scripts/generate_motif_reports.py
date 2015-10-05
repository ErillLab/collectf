# This script generates motif reports for all TFs, species and experimental
# technique levels. The generated motif reports are saved and used for browse
# pages in CollecTF without accessing database and computing motif reports for
# the query TFs, species or experimental techniques.

import os

from tqdm import tqdm
import pickle

from base import models
from browse import motif_report
from collectf import settings

def pickle_report(reports, by, object_id):
    """Pickles the given report.

    The pickling type can be taxonomy, TF, TF_family, experimental technique and
    experimental technique category.
    """
    assert by in ['taxonomy', 'TF', 'TF_family', 'experimental_technique',
                  'experimental_technique_category']
    pickle_file = os.path.join(settings.PICKLE_ROOT, 'motif_reports',
                               by + '_' + str(object_id) + '.pkl')
    pickle.dump(reports, open(pickle_file, 'w'))

def motif_reports_by_tf_family():
    """Generates motif reports for each TF-family."""
    TF_families = models.TFFamily.objects.all()
    for TF_family in tqdm(TF_families):
        TFs = models.TF.objects.filter(family=TF_family)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF__in=TFs)
        reports = motif_report.make_reports(curation_site_instances)
        pickle_report(reports, 'TF_family', TF_family.TF_family_id)

def motif_reports_by_tf():
    """Generates motif reports for each TF."""
    for TF in tqdm(models.TF.objects.all()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF=TF)
        reports = motif_report.make_reports(curation_site_instances)
        pickle_report(reports, 'TF', TF.TF_id)

def motif_reports_by_taxonomy():
    """Generates motif reports for each taxon."""
    for taxon in tqdm(models.Taxonomy.objects.all()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            site_instance__genome__taxonomy__in=taxon.get_all_species())
        reports = motif_report.make_reports(curation_site_instances)
        pickle_report(reports, 'taxonomy', taxon.taxonomy_id)

def motif_reports_by_experimental_technique_category():
    """Generate motif reports for experimental technique categories."""
    for category in tqdm(models.ExperimentalTechniqueCategory.objects.all()):
        techniques = models.ExperimentalTechnique.objects.filter(
            categories=category)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            experimental_techniques__in=techniques)
        reports = motif_report.make_reports(curation_site_instances)
        pickle_report(reports, 'experimental_technique_category',
                      category.category_id)
    
def motif_reports_by_experimental_technique():
    """Generates motif reports for each experimental technique."""
    for technique in tqdm(models.ExperimentalTechnique.objects.all()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            experimental_techniques=technique)
        reports = motif_report.make_reports(curation_site_instances)
        pickle_report(reports, 'experimental_technique', technique.technique_id)

def run():
    print "Generating motif reports for all experimental techniques."
    motif_reports_by_experimental_technique()
    
    print "Generating motif reports for all experimental technique categories."
    motif_reports_by_experimental_technique_category()

    print "Generating motif reports for all taxonomy levels."""
    motif_reports_by_taxonomy()
    
    print "Generating motif reports for TF families."
    motif_reports_by_tf_family()

    print "Generating motif reports for TFs."
    motif_reports_by_tf()


