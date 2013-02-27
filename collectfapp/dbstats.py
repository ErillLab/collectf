"""Set of functions for database statistics"""

from models import *

def curation_stats(csv_file="db_stats.csv"):
    tab = '\t'
    newline = '\n'
    header = tab.join(("Curation id",
                       "Publication id",
                       "Curator",
                       "TF instance",
                       "TF",
                       "TF family",
                       "TF function",
                       "TF type",
                       "TF species",
                       "site species",
                       "annotated_sites",
                       "sites_genome_accession",
                       ))

    f = open(csv_file, 'w')
    f.write(header + newline)
    for c in Curation.objects.all():
        curation_id = c.curation_id
        publication_id = c.publication.publication_id
        curator = c.curator.user.username
        TF_instance = c.TF_instance.name
        TF = c.TF.name
        TF_family = c.TF.family.name
        TF_function = c.TF_function
        TF_type = c.TF_type
        TF_species = c.TF_species
        site_species = c.site_species

        # all site instances should belong to the same genome!
        site_instances = c.site_instances.all()
        assert len(site_instances) > 0, "no site instance"
        assert all(x.genome.genome_accession ==
                   site_instances[0].genome.genome_accession
                   for x in site_instances)

        annotated_sites = len(c.site_instances.all())
        sites_genome_accession = c.site_instances.all()[0].genome.genome_accession

        # write to csv
        f.write(tab.join(map(str, (curation_id,
                                   publication_id,
                                   curator,
                                   TF_instance,
                                   TF,
                                   TF_family,
                                   TF_function,
                                   TF_type,
                                   TF_species,
                                   site_species,
                                   annotated_sites,
                                   sites_genome_accession))))
        f.write(newline)
    f.close()
            

    
def exp_technique_stats():
    """Frequency of each experimental techniques used in curations"""
    all_curations = Curation.objects.all()
    exp_techniques_freq = dict((e.name, 0) for e in ExperimentalTechnique.objects.all())
    for c in all_curations:
        for e in c.experimental_techniques.all():
            exp_techniques_freq[e.name] += 1

    for k,v in exp_techniques_freq.items():
        print k,v
        
