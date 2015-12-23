"""Generate GO annotations."""

import os
import time

from django.conf import settings

from core import models


def new_gpad_record(uniprot_id, relationship, go_id, db_ref, evidence_code):
    return dict(db='UniProtKB',
                db_object_id=uniprot_id,
                relationship=relationship,
                go_id=go_id,
                db_ref=('PMID:' + db_ref),
                evidence_code=evidence_code,
                with_from='',
                interacting_taxon_id='',
                date=time.strftime('%Y%m%d'),
                assigned_by='CollecTF',
                annotation_extension='',
                annotation_properties='')


def gpad_record_to_str(gpad_record):
    return '\t'.join(gpad_record[field] for field in
                     ['db', 'db_object_id', 'relationship', 'go_id', 'db_ref',
                      'evidence_code', 'with_from', 'interacting_taxon_id',
                      'date', 'assigned_by', 'annotation_extension',
                      'annotation_properties'])


def go_annotations(TF_instance):
    """Generates list of GO annotations for the given TF instance."""
    def get_techniques(regulations):
        """Gets list of techniques for a given set of regulations."""
        technique_ids = regulations.values_list(
            'curation_site_instance__experimental_techniques',
            flat=True).distinct()
        return models.ExperimentalTechnique.objects.filter(
            pk__in=technique_ids, EO_term__isnull=False)

    print TF_instance
    go_records = []

    # TF-binding-annotations (enables) molecular function
    regulations = models.Regulation.objects.filter(
        evidence_type='exp_verified',
        curation_site_instance__curation__TF_instances=TF_instance)
    if regulations.filter(
            curation_site_instance__TF_function__in=['ACT', 'REP']):
        # If there is evidence of regulation AND at least one gene is shown to
        # be activated or repressed,
        for regulation in regulations.filter(
                curation_site_instance__TF_function='ACT'):
            for tech in regulation.binding_experimental_techniques:
                if not tech.EO_term:
                    continue
                go_records.append(new_gpad_record(
                    uniprot_id=TF_instance.uniprot_accession,
                    relationship='enables', go_id='GO:0001216',
                    db_ref=regulation.ref_pmid, evidence_code=tech.EO_term))
        for regulation in regulations.filter(
                curation_site_instance__TF_function='REP'):
            for tech in regulation.binding_experimental_techniques:
                if not tech.EO_term:
                    continue
                go_records.append(new_gpad_record(
                    uniprot_id=TF_instance.uniprot_accession,
                    relationship='enables', go_id='GO:0001217',
                    db_ref=regulation.ref_pmid, evidence_code=tech.EO_term))
        for regulation in regulations:
            for tech in regulation.binding_experimental_techniques:
                if not tech.EO_term:
                    continue
                go_records.append(new_gpad_record(
                    uniprot_id=TF_instance.uniprot_accession,
                    relationship='enables', go_id='GO:0000976',
                    db_ref=regulation.ref_pmid, evidence_code=tech.EO_term))
    elif regulations:
        # If there is evidence of TF-mediated regulation for at least one gene,
        # but no indication of whether the TF up- or down-regulates
        # (repressor/activator) for any gene in the curation
        pass

    return set(map(gpad_record_to_str, go_records))


def generate_gpad_file():
    """Collects GO annotations for all TF-instances and generates GPAD file."""
    header = [
        "!gpa-version: 1.1",
        "!Submission Date: 12/22/2015",
        "!",
        "!Project_name: CollecTF - Experimentally validated",
        "!URL: http://www.collectf.org/",
        "!Contact Email: collectf@umbc.edu",
        "!",
        '\t'.join(['!DB', 'DB_Object_ID', 'Relationship', 'GO ID', 'Reference',
                   'Evidence code', 'With (or) From', 'Interacting taxon ID',
                   'Date', 'Assigned By'])]
    content = [line
               for TF_instance in models.TFInstance.objects.all()
               for line in go_annotations(TF_instance)]
    export_file = os.path.join(settings.STATICFILES_DIRS[0], 'collectf.gpad')
    with open(export_file, 'w') as f:
        f.write('\n'.join(header + content))


def run():
    generate_gpad_file()
