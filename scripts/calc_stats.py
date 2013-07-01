import pickle
from collectfapp import models

def curation_stats():
    """Curation statistics.. Count the number of curations/sites for each TF and
    species in the database.
    """
    all_TFs = models.TF.objects.all()
    all_species = models.Strain.objects.all()

    num_curations_by_TF_species = {}
    num_sites_by_TF_species = {}
    
    for TF in all_TFs:
        num_curations_by_TF_species[TF.name] = {}
        num_sites_by_TF_species[TF.name] = {}
        for sp in all_species:
            # get curation_site_instance objects
            csi = models.Curation_SiteInstance.objects.filter(
                site_instance__genome__strain=sp,
                curation__TF=TF)
            num_sites = csi.values_list('site_instance', flat=True).distinct().count()
            num_curations = csi.values_list('curation', flat=True).distinct().count()
            num_curations_by_TF_species[TF.name][sp.name] = num_curations
            num_sites_by_TF_species[TF.name][sp.name] = num_sites

    master_dict = dict(
        num_TFs = len(all_TFs),
        num_species = len(all_species),
        num_curations = len(models.Curation.objects.all()),
        num_sites = len(models.SiteInstance.objects.all()),
        num_publications = len(models.Publication.objects.all()),
        pub_completed = '%.1f' % (len(models.Publication.objects.filter(curation_complete=True)) * 100.0 /
                                  len(models.Publication.objects.all())),
        num_curations_by_TF_species = num_curations_by_TF_species,
        num_sites_by_TF_species = num_sites_by_TF_species,
        TFs= [tf.name for tf in all_TFs],
        species = [sp.name for sp in all_species]
        )

    # write stats to pickle
    pickle.dump(master_dict, open('dbstats.pickle', 'w'))

def run():
    curation_stats()
    print "I am a script!"
