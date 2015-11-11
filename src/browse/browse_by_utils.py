def curation_site_instances_values_list(curation_site_instances):
    """Returns the values list of Curation_SiteInstance objects.

    https://docs.djangoproject.com/en/1.8/ref/models/querysets/#values-list
    """
    return curation_site_instances.values_list(
        'curation__TF_instances__TF__TF_id',
        'curation__TF_instances__TF__name',
        'site_instance__genome__taxonomy__pk',
        'site_instance__genome__taxonomy__name').distinct().order_by(
            'curation__TF_instances__TF__name',
            'site_instance__genome__taxonomy__name')
