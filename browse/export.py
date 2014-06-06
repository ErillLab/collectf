"""Export functions to export binding site data in various formats."""
import models
from base import bioutils
from django.http import HttpResponse

def export_sites(request):
    """Given a list of sites, report fasta/csv/arff etc. file containing sites
    for particular TF and species"""

    export_options = ['fasta', 'tsv', 'tsv-raw', 'arff',
                      'PSFM-jaspar', 'PSFM-transfac', 'PSFM-raw-fasta']
    ext = {
        'fasta': 'fas',
        'tsv': 'tsv',
        'tsv-raw': 'tsv',
        'arff': 'arff',
        'PSFM-jaspar': 'mat',
        'PSFM-transfac': 'mat',
        'PSFM-raw-fasta': 'mat'
    }
    func = {
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

    # Given the list of curation-site-instance objects,
    curation_site_id_groups = request.POST.getlist('site_id')
    meta_sites = [models.Curation_SiteInstance.objects.filter(pk__in=id_group.split('|'))
                  for id_group in curation_site_id_groups]

    response = HttpResponse(content_type='application/download')
    # call corresponding export function
    response["Content-Disposition"] = \
            "attachment;filename=collectf-export-%s.%s" % (export_format, ext[export_format])
    response.write(func[export_format](meta_sites, **kwargs[export_format]))
    return response

def export_base(meta_sites):
    """Base function for export. Called by format-specific export functions."""
    rows = []
    for i,meta_site in enumerate(meta_sites):
        values = meta_site.values('curation__TF__name',
                                  'curation__TF_instances__protein_accession',
                                  'site_instance__genome__genome_accession',
                                  'site_instance__genome__organism').distinct()
        assert len(values)==1, values
        values = values[0]
        values['start_pos'] = meta_site[0].site_instance.start+1
        values['end_pos'] = meta_site[0].site_instance.end+1
        values['strand'] = meta_site[0].site_instance.strand
        values['seq'] = meta_site[0].site_instance.seq
        values['mode'] = meta_site[0].TF_function
        rows.append(values)
    return rows

def export_fasta(meta_sites, **kwargs):
    """Export TFBS in FASTA format."""
    fasta_str = ""
    # dont load genome sequence when retrieving from DB
    rows = export_base(meta_sites)
    for row in rows:
        fasta_str += (">TF_%(TF)s_%(TF_accession)s|genome_%(genome)s_%(organism)s|start=%(start_pos)d|end=%(end_pos)d|strand=%(strand)d\n" %
                      {'TF': row['curation__TF__name'],
                       'TF_accession': row['curation__TF_instances__protein_accession'],
                       'genome': row['site_instance__genome__genome_accession'],
                       'organism': row['site_instance__genome__organism'],
                       'start_pos': row['start_pos'],
                       'end_pos': row['end_pos'],
                       'strand': row['strand'],
                       })
        fasta_str += ("%s\n" % row['seq'])
    return fasta_str

tsv_sep = '\t'
tsv_header = tsv_sep.join(['TF',
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
    """Export tab-separated file."""
    tsv_lines = []
    sep = '\t'
    tsv_lines.append(tsv_header)
    rows = export_base(meta_sites)
    for i,row in enumerate(rows):
        csis = meta_sites[i]
        experimental_evidence = [csi.experimental_techniques.all() for csi in csis.all()]
        pmids = [csi.curation.publication.pmid for csi in csis.all()]
        regulated_genes = models.Regulation.objects.filter(evidence_type="exp_verified")\
                          .filter(curation_site_instance__in=csis)
        line = tsv_sep.join(map(str, [row['curation__TF__name'],
                                      row['curation__TF_instances__protein_accession'],
                                      row['site_instance__genome__genome_accession'],
                                      row['site_instance__genome__organism'],
                                      row['start_pos'],
                                      row['end_pos'],
                                      row['strand'],
                                      row['seq'],
                                      row['mode'],
                                      ' | '.join((','.join(t.name for t in evidence) + ' [PMID:%s]' % pmid)
                                                 for evidence,pmid in zip(experimental_evidence, pmids)),
                                      ', '.join(reg.gene.locus_tag for reg in regulated_genes),
                                      ]))
        tsv_lines.append(line)
    return '\n'.join(tsv_lines)

def export_tsv_raw(meta_sites, **kwargs):
    tsv_lines = []
    tsv_lines.append(tsv_header)
    rows = export_base(meta_sites)
    for i,row in enumerate(rows):
        csis = meta_sites[i]
        # report each csi individually
        for csi in csis.all():
            line = tsv_sep.join(map(str, [csi.curation.TF.name,
                                          csi.curation.TF_instances.all()[0].protein_accession,
                                          csi.site_instance.genome.genome_accession,
                                          csi.site_instance.genome.organism,
                                          csi.site_instance.start,
                                          csi.site_instance.end,
                                          csi.site_instance.strand,
                                          csi.site_instance.seq,
                                          csi.TF_function,
                                          ', '.join(t.name for t in csi.experimental_techniques.all()) + \
                                                                       ' [PMID:%s]' % csi.curation.publication.pmid,
                                          ','.join(r.gene.locus_tag for r in models.Regulation.objects.filter(evidence_type="exp_verified")\
                                                                                                      .filter(curation_site_instance=csi)),
                                          ]))
            tsv_lines.append(line)
    return '\n'.join(tsv_lines)


def export_arff(meta_sites, **kwargs):
    arff_lines = []
    # comments
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
    for i,row in enumerate(rows):
        arff_lines.append(','.join([row['curation__TF__name'],
                                    row['curation__TF_instances__protein_accession'],
                                    row['site_instance__genome__genome_accession'],
                                    '"' + row['site_instance__genome__organism'] + '"',
                                    str(row['start_pos']),
                                    str(row['end_pos']),
                                    str(row['strand']),
                                    row['seq']]))
    return '\n'.join(arff_lines)

def export_PSFM(meta_sites, **kwargs):
    """Export Position-Specific-Frequency-Matrix"""
    format = kwargs['format']
    rows = export_base(meta_sites)
    aligned = bioutils.run_lasagna([m[0].site_instance for m in meta_sites])
    motif = bioutils.build_motif(aligned)
    consensus = bioutils.degenerate_consensus(motif)

    TF_name= ','.join(set(row['curation__TF__name'] for row in rows))
    sp = ','.join(set('_'.join(row['site_instance__genome__organism'].split()) for row in rows))
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
