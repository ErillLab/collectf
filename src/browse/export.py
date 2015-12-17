"""Export functions to export binding site data in various formats."""

from django.http import HttpResponse
from django.views.decorators.http import require_POST

from core import bioutils
from core import models
from core import metasite
from core import lasagna

@require_POST
def export_sites(request):
    """Exports report in different formats for the given list of sites."""
    export_options = ['fasta', 'tsv', 'tsv-raw', 'arff', 'PSFM-jaspar',
                      'PSFM-transfac', 'PSFM-raw-fasta']
    extensions = {
        'fasta': 'fas',
        'tsv': 'tsv',
        'tsv-raw': 'tsv',
        'arff': 'arff',
        'PSFM-jaspar': 'mat',
        'PSFM-transfac': 'mat',
        'PSFM-raw-fasta': 'mat'
    }
    export_functions = {
        'fasta': export_fasta,
        'tsv': export_tsv,
        'tsv-raw': export_tsv_raw,
        'arff': export_arff,
        'PSFM-jaspar': export_PSFM,
        'PSFM-transfac': export_PSFM,
        'PSFM-raw-fasta': export_PSFM,
    }
    kwargs = {
        'fasta': {},
        'tsv': {},
        'tsv-raw': {},
        'arff': {},
        'PSFM-jaspar': {'format':'JASPAR'},
        'PSFM-transfac': {'format':'TRANSFAC'},
        'PSFM-raw-fasta': {'format': 'raw_fasta'},
    }

    for opt in export_options:
        if opt in request.POST:
            export_format = opt
            break
    else: assert False, "invalid export option"

    meta_sites = metasite.create_meta_sites(
        models.Curation_SiteInstance.objects.filter(
            pk__in=request.POST.getlist('site_id')))
                                   
    response = HttpResponse(content_type='application/download')

    # Call corresponding export function
    response["Content-Disposition"] = (
        "attachment;filename=collectf-export-%s.%s" %
        (export_format, extensions[export_format]))
    response.write(export_functions[export_format](
        meta_sites, **kwargs[export_format]))
    return response

def export_base(meta_sites):
    """Helper function to be called by export functions.

    Called by format-specific export functions.
    """
    rows = []
    for meta_site in meta_sites:
        delegate = meta_site.delegate
        rows.append({
            'TF_name': delegate.curation.TF.name,
            'TF_accession': delegate.curation.TF_instance_accessions[0],
            'genome_accession': delegate.site_instance.genome.genome_accession,
            'organism': delegate.site_instance.genome.organism,
            'start_pos': delegate.site_instance.start+1,
            'end_pos': delegate.site_instance.end+1,
            'strand': delegate.site_instance.strand,
            'seq': delegate.site_instance.sequence,
            'mode': delegate.TF_function
        })
    return rows

def export_fasta(meta_sites, **kwargs):
    """Exports binding sites in FASTA format."""
    fasta_str = ""
    rows = export_base(meta_sites)
    for row in rows:
        fasta_str += (
            '>TF_%(TF)s_%(TF_accession)s|genome_%(genome)s %(organism)s|'
            'start=%(start_pos)d|end=%(end_pos)d|strand=%(strand)d\n' %
            {'TF': row['TF_name'],
             'TF_accession': row['TF_accession'],
             'genome': row['genome_accession'],
             'organism': row['organism'],
             'start_pos': row['start_pos'],
             'end_pos': row['end_pos'],
             'strand': row['strand']})
        fasta_str += ('%s\n' % row['seq'])
    return fasta_str

TSV_SEPARATOR = '\t'
TSV_HEADER = TSV_SEPARATOR.join(['TF',
                                 'TF_accession',
                                 'genome_accession',
                                 'organism',
                                 'site_start',
                                 'site_end',
                                 'site_strand',
                                 'sequence',
                                 'mode',
                                 'experimental_evidence',
                                 'regulated genes (locus_tags)'
])

def export_tsv(meta_sites, **kwargs):
    """Exports tab-separated file."""
    tsv_lines = []
    tsv_lines.append(TSV_HEADER)
    rows = export_base(meta_sites)
    for meta_site, row in zip(meta_sites, rows):
        experimental_evidence = [csi.experimental_techniques.all()
                                 for csi in meta_site.curation_site_instances]
        pmids = [csi.curation.publication.pmid
                 for csi in meta_site.curation_site_instances]
        regulated_genes = models.Regulation.objects.filter(
            evidence_type='exp_verified').filter(
                curation_site_instance__in=meta_site.curation_site_instances)
        line = TSV_SEPARATOR.join([
            row['TF_name'],
            row['TF_accession'],
            row['genome_accession'],
            row['organism'],
            str(row['start_pos']),
            str(row['end_pos']),
            str(row['strand']),
            row['seq'],
            row['mode'],
            ' | '.join(
                (', '.join(t.name for t in evidence) + ' [PMID:%s]' % pmid)
                for evidence, pmid in zip(experimental_evidence, pmids)),
            ', '.join(reg.gene.locus_tag for reg in regulated_genes)])
        tsv_lines.append(line)
    return '\n'.join(tsv_lines)

def export_tsv_raw(meta_sites, **kwargs):
    """Exports each curation-site-instance individually."""
    tsv_lines = []
    tsv_lines.append(TSV_HEADER)
    for meta_site in meta_sites:
        # Report each csi individually
        for csi in meta_site.curation_site_instances:
            line = TSV_SEPARATOR.join([
                csi.curation.TF.name,
                csi.curation.TF_instance_accessions[0],
                csi.site_instance.genome.genome_accession,
                csi.site_instance.genome.organism,
                str(csi.site_instance.start+1),
                str(csi.site_instance.end+1),
                str(csi.site_instance.strand),
                csi.site_instance.sequence,
                csi.TF_function,
                (', '.join(t.name for t in csi.experimental_techniques.all()) +
                 ' [PMID:%s]' % csi.curation.publication.pmid),
                ', '.join(r.gene.locus_tag
                          for r in models.Regulation.objects.filter(
                                  evidence_type="exp_verified",
                                  curation_site_instance=csi))])
            tsv_lines.append(line)
    return '\n'.join(tsv_lines)


def export_arff(meta_sites, **kwargs):
    """Exports sites in arff format -- mainly used by Weka."""
    arff_lines = []
    # Comments
    arff_lines.append('@RELATION "CollecTF export"')
    arff_lines.append('@ATTRIBUTE TF_name STRING')
    arff_lines.append('@ATTRIBUTE Protein_accession STRING')
    arff_lines.append('@ATTRIBUTE Genome_accession STRING')
    arff_lines.append('@ATTRIBUTE Genome_organism STRING')
    arff_lines.append('@ATTRIBUTE TF_binding_site_start NUMERIC')
    arff_lines.append('@ATTRIBUTE TF_binding_site_end NUMERIC')
    arff_lines.append('@ATTRIBUTE TF_binding_site_strand NUMERIC')
    arff_lines.append('@ATTRIBUTE TF_binding_site_sequence STRING')

    arff_lines.append('@DATA')

    rows = export_base(meta_sites)
    for row in rows:
        arff_lines.append(','.join([
            row['TF_name'],
            row['TF_accession'],
            row['genome_accession'],
            '"' + row['organism'] + '"',
            str(row['start_pos']),
            str(row['end_pos']),
            str(row['strand']),
            row['seq']]))
    return '\n'.join(arff_lines)


def export_PSFM(meta_sites, **kwargs):
    """Exports Position-Specific-Frequency-Matrix"""
    format = kwargs['format']
    rows = export_base(meta_sites)
    aligned = lasagna.lasagna([m.delegate_site_instance for m in meta_sites])
    motif = bioutils.build_motif(aligned)
    consensus = motif.degenerate_consensus

    TF_name= ','.join(set(row['TF_name'] for row in rows))
    sp = ','.join(set('_'.join(row['organism'].split()) for row in rows))
    lines = []
    if format == 'JASPAR':
        lines.append('> CollecTF_%s_%s' % (TF_name, sp))
        lines.append('A [ %s ]' % (' '.join(map(str, motif.counts['A']))))
        lines.append('C [ %s ]' % (' '.join(map(str, motif.counts['C']))))
        lines.append('G [ %s ]' % (' '.join(map(str, motif.counts['G']))))
        lines.append('T [ %s ]' % (' '.join(map(str, motif.counts['T']))))
    elif format == 'TRANSFAC':
        lines.append('ID %s' % TF_name)
        lines.append('BF %s' % sp)
        lines.append('PO\tA\tC\tG\tT')
        lines.extend('%02d\t%d\t%d\t%d\t%d\t%s' % (po+1, motif.counts['A'][po],
                                                   motif.counts['C'][po],
                                                   motif.counts['G'][po],
                                                   motif.counts['T'][po],
                                                   consensus[po])
                     for po in range(motif.length))
        lines.append('XX')
    elif format == 'raw_fasta':
        lines.append('>CollecTF_%s_%s' % (TF_name, sp))
        lines.extend('%d\t%d\t%d\t%d' % (motif.counts['A'][po],
                                         motif.counts['C'][po],
                                         motif.counts['G'][po],
                                         motif.counts['T'][po])
                     for po in range(motif.length))

    return '\n'.join(lines)

