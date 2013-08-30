from collectfapp import models
from collectfapp import sitesearch

def fix_operon():
    for csi in models.Curation_SiteInstance.objects.iterator():
        regs = csi.regulation_set.all()
        site = csi.site_instance
        genes = models.Gene.objects.filter(genome=site.genome).order_by('start')
        nearby_genes = sitesearch.locate_nearby_genes(genes, site)

        correct_nearby = []
        incorrect_nearby = []
        for reg in regs:
            if reg.gene in nearby_genes:
                correct_nearby.append(reg)
            else:
                incorrect_nearby.append(reg)
                if reg.evidence_type!='inferred':
                    print csi

def run():
    fix_operon()
    print 'fixed operons'
