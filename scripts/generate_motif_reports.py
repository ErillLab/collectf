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

def pickle_report(curation_site_instances, filename_suffix):
    """Generates a report and pickles it, given Curation_SiteInstances objects.
    """
    reports = motif_report.make_reports(curation_site_instances)
    reports_file = os.path.join(
        settings.PICKLE_ROOT, 'reports', filename_suffix + '.pkl')
    pickle.dump(reports, open(reports_file, 'w'))
    
    ensemble_report = motif_report.make_ensemble_report(curation_site_instances)
    ensemble_report_file = os.path.join(
        settings.PICKLE_ROOT, 'ensemble_reports', filename_suffix + '.pkl')
    pickle.dump(ensemble_report, open(ensemble_report_file, 'w'))

def motif_reports_by_tf_family():
    """Generates motif reports for each TF-family."""
    TF_families = models.TFFamily.objects.all()
    for TF_family in tqdm(TF_families):
        TFs = models.TF.objects.filter(family=TF_family)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF__in=TFs)
        pickle_report(curation_site_instances,
                      'TF_family_%s' % TF_family.TF_family_id)

def motif_reports_by_tf():
    """Generates motif reports for each TF."""
    for TF in tqdm(models.TF.objects.all()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF=TF)
        pickle_report(curation_site_instances, 'TF_%s' % TF.TF_id)

def motif_reports_by_taxonomy():
    """Generates motif reports for each taxon."""
    for taxon in tqdm(models.Taxonomy.objects.all()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            site_instance__genome__taxonomy__in=taxon.get_all_species())
        pickle_report(curation_site_instances,
                      'taxonomy_%s' % taxon.taxonomy_id)

def motif_reports_by_technique_function():
    """Generates motif reports for all binding/expression techniques."""
    for function in ['binding', 'expression']:
        techniques = models.ExperimentalTechnique.objects.filter(
            preset_function=function)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            experimental_techniques__in=techniques)
        pickle_report(curation_site_instances,
                      'experimental_technique_all_' + function)

def motif_reports_by_experimental_technique_category():
    """Generate motif reports for experimental technique categories."""
    # TODO(sefa): Group motifs by binding/expression function as well.
    for category in tqdm(models.ExperimentalTechniqueCategory.objects.all()):
        techniques = models.ExperimentalTechnique.objects.filter(
            categories=category)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            experimental_techniques__in=techniques)
        pickle_report(
            curation_site_instances,
            'experimental_technique_category_%s' % category.category_id)
    
def motif_reports_by_experimental_technique():
    """Generates motif reports for each experimental technique."""
    for technique in tqdm(models.ExperimentalTechnique.objects.all()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            experimental_techniques=technique)
        pickle_report(curation_site_instances,
                      'experimental_technique_%s' % technique.technique_id)

def motif_reports_by_species_and_TF():
    """Generates motif reports for each pair of species and TF."""
    TF_and_species = models.Curation_SiteInstance.objects.values_list(
        'curation__TF_instances__TF',
        'site_instance__genome__taxonomy')
    for TF, species in tqdm(TF_and_species.distinct()):
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF=TF,
            site_instance__genome__taxonomy=species)
        pickle_report(curation_site_instances,
                      'tf_%d_species_%d' % (TF, species))

def run():
    print "Generating motif reports for all pairs of TFs and species."
    motif_reports_by_species_and_TF()
    
    print "Generating motif reports for all experimental techniques."
    motif_reports_by_experimental_technique()
    
    print "Generating motif reports for all experimental technique categories."
    motif_reports_by_experimental_technique_category()

    print "Generating motif reports for all binding/expression techniques."
    motif_reports_by_technique_function()

    print "Generating motif reports for all taxonomy levels."""
    motif_reports_by_taxonomy()
    
    print "Generating motif reports for TF families."
    motif_reports_by_tf_family()

    print "Generating motif reports for TFs."
    motif_reports_by_tf()




