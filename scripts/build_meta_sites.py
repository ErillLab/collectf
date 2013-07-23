import pickle
import os
from collectfapp import models
from collectf import settings

def get_overlap(a, b):
    return max(0, min(a[1], b[1]) -  max(a[0], b[0]))

def get_meta_site_loc(ms):
    # given a meta site find its start and end positions
    start = min(map(lambda s: s.site_instance.start, ms.curation_siteinstance_set.all()))
    end = max(map(lambda s: s.site_instance.end, ms.curation_siteinstance_set.all()))
    return start,end

def build_meta_sites():
    # get all site instances
    for csi in models.Curation_SiteInstance.objects.filter(is_motif_associated=True).iterator():
        loc_a = (csi.site_instance.start, csi.site_instance.end)
        genome = csi.site_instance.genome
        TF_instance = csi.curation.TF_instance
        for ms in models.MetaSiteInstance.objects.filter(genome=genome, TF_instance=TF_instance).iterator():
            loc_b = get_meta_site_loc(ms)
            # if meta-site overlaps site-instance >75% of site-len, add
            if get_overlap(loc_a, loc_b) >= 0.75 * (loc_a[1]-loc_a[0]):
                csi.meta_site_instance = ms
                csi.save()
                ms.start = min(loc_a[0], loc_b[0])
                ms.end = max(loc_a[1], loc_b[1])
                ms.save()
                print 'found'
                break
        else: # in case meta-site is not found
            ms = models.MetaSiteInstance(genome=genome,
                                         TF_instance=TF_instance,
                                         start=csi.site_instance.start,
                                         end=csi.site_instance.end)
            ms.save()
            csi.meta_site_instance = ms
            csi.save()
            print 'created new'

def find_meta_sites(with_at_least=2):
    for ms in models.MetaSiteInstance.objects.iterator():
        if ms.curation_siteinstance_set.values('site_instance').distinct().count() >= with_at_least:
            print ms.pk, ms.curation_siteinstance_set.all()

def run():
    #build_meta_sites()
    find_meta_sites()
    print "I am a script!"
