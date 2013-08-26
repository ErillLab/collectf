from browse_base import *

def export_sites(request):
    """Given a list of sites, report FASTA/CSV file
    containing sites for particular TF and species"""

    if 'fasta' in request.POST: export_format = 'fasta'
    elif 'tsv'in request.POST: export_format = 'tsv'
    elif 'tsv-raw' in request.POST: export_format='tsv_raw'
    assert export_format

    # given the list of curation-site-instance objects,
    curation_site_id_groups = request.POST.getlist('site_id')

    meta_sites = [models.Curation_SiteInstance.objects.filter(pk__in=id_group.split('|')) for id_group in curation_site_id_groups]
    
    response = HttpResponse(content_type='application/download')
    if export_format=="fasta":
        response['Content-Disposition'] = 'attachment;filename=collectf_export.fas'
        response.write(export_fasta(meta_sites))
    elif export_format=="tsv":
        response['Content-Disposition'] = 'attachment;filename=collectf_export.tsv'
        response.write(export_tsv(meta_sites))
    elif export_format=="tsv_raw":
        response['Content-Disposition'] = 'attachment;filename=collectf_export.tsv'
        response.write(export_tsv_raw(meta_sites))
    return response



def export_base(meta_sites):
    rows = []
    for meta_site in meta_sites:
        values = meta_site.values('curation__TF__name',
                                  'curation__TF_instance__protein_accession',
                                  'site_instance__genome__genome_accession',
                                  'site_instance__genome__organism').distinct()
        
        assert len(values)==1
        values = values[0]
        values['start_pos'] = meta_site[0].site_instance.start+1
        values['end_pos'] = meta_site[0].site_instance.end+1
        values['strand'] = meta_site[0].site_instance.strand
        values['seq'] = meta_site[0].site_instance.seq
        rows.append(values)
    return rows
    
    
def export_fasta(meta_sites):
    fasta_str = ""
    # dont load genome sequence when retrieving from DB
    rows = export_base(meta_sites)
    for row in rows:
        fasta_str += (">TF_%(TF)s_%(TF_accession)s|genome_%(genome)s_%(organism)s|start=%(start_pos)d|end=%(end_pos)d|strand=%(strand)d\n" %
                      {'TF': row['curation__TF__name'],
                       'TF_accession': row['curation__TF_instance__protein_accession'],
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
                           'experimental_evidence',
                           'regulated genes (locus_tags)'
                           ])

def export_tsv(meta_sites):
    tsv_lines = []
    sep = '\t'
    tsv_lines.append(tsv_header)
    rows = export_base(meta_sites)
    for i,row in enumerate(rows):
        csis = meta_sites[i]
        experimental_evidence = [csi.curation.experimental_techniques.all() for csi in csis.all()]
        pmids = [csi.curation.publication.pmid for csi in csis.all()]
        regulated_genes = models.Regulation.objects.filter(evidence_type="exp_verified")\
                          .filter(curation_site_instance__in=csis)
        line = tsv_sep.join(map(str, [row['curation__TF__name'],
                                 row['curation__TF_instance__protein_accession'],
                                 row['site_instance__genome__genome_accession'],
                                 row['site_instance__genome__organism'],
                                 row['start_pos'],
                                 row['end_pos'],
                                 row['strand'],
                                 row['seq'],
                                 ' | '.join((','.join(t.name for t in evidence) + '[PMID:%s]' % pmid) for evidence,pmid in zip(experimental_evidence, pmids)),
                                 ','.join(reg.gene.locus_tag for reg in regulated_genes),
                                 ]))
        tsv_lines.append(line)
    return '\n'.join(tsv_lines)

def export_tsv_raw(meta_sites):
    tsv_lines = []
    tsv_lines.append(tsv_header)
    rows = export_base(meta_sites)
    for i,row in enumerate(rows):
        csis = meta_sites[i]
        # report each csi individually
        for csi in csis.all():
            line = tsv_sep.join(map(str, [csi.curation.TF.name,
                                          csi.curation.TF_instance.protein_accession,
                                          csi.site_instance.genome.genome_accession,
                                          csi.site_instance.genome.organism,
                                          csi.site_instance.start,
                                          csi.site_instance.end,
                                          csi.site_instance.strand,
                                          csi.site_instance.seq,
                                          ','.join(t.name for t in csi.curation.experimental_techniques.all()) + '[PMID:%s]' % csi.curation.publication.pmid,
                                          ','.join(r.gene.locus_tag for r in models.Regulation.objects.filter(evidence_type="exp_verified")\
                                                                                                      .filter(curation_site_instance=csi)),
                                          ]))
            tsv_lines.append(line)
    return '\n'.join(tsv_lines)
                                    




