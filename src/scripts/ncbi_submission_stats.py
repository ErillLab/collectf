# NCBI submission stats

from core import models


def ncbi_submission_stats():
    """Prints the NCBI submission statistics."""

    ncbi_submissions = models.NCBISubmission.objects
    papers = ncbi_submissions.values_list(
        'curation_site_instance__curation__publication').distinct()
    print "Annotated papers:", papers.count()
    TFs = ncbi_submissions.values_list(
        'curation_site_instance__curation__TF_instances_TF').distinct()
    print "Annotated TFs:", TFs.count()
    sites = ncbi_submissions.values_list(
        'curation_site_instance__site_instance').distinct()
    print "Annotated sites:", sites.count()
    genomes = ncbi_submissions.values_list(
        'curation_site_instance__site_instance__genome').distinct()
    print "Annotated genomes:", genomes.count()


def run():
    ncbi_submission_stats()
