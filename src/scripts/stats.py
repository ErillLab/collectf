from core import models

def TF_stats():
    #TF_instance_ids = models.Curation.objects.values_list('TF_instances').distinct()
    TF_instance_ids = models.Curation_SiteInstance.objects.values_list('curation__TF_instances').distinct()
    

    for TF in models.TF.objects.all():
        #print TF.name
        print models.TFInstance.objects.filter(
            TF_instance_id__in=TF_instance_ids,
           TF=TF).count()

            

def run():
    TF_stats()
        
        
