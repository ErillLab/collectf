import models
import re
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.shortcuts import render_to_response
from django.template import RequestContext
import collectfapp
import views
import utils

@login_required
def export_db(request):
    return export2tsv(request)
    return render_to_response("export_db.html",
                              context_instance=RequestContext(request))
    

def export2tsv(request):
    filename = 'collectf_db_export.tab'
    separator = '\t'
    tsv_file = HttpResponse(content_type='application/download')
    tsv_file['Content-Disposition'] = 'attachment;filename="%s"' % filename

    tsv_file.write(separator.join(['genome',
                                   'TF',
                                   'start_pos',
                                   'end_pos',
                                   'strand',
                                   'sequence',
                                   'regulated_genes (locus_tag)',
                                   'PMID',
                                   'Experimental evidence',
                                   'DBXRef',
                                   'DBXRef_link'
                                   ]))
    tsv_file.write('\n')

    for i,csi in enumerate(models.Curation_SiteInstance.objects.filter(curation__NCBI_submission_ready=True).iterator()):
        print i
        fields = []
        fields.append(csi.site_instance.genome.genome_accession) # genome accession
        fields.append(csi.curation.TF.name +
                      ' [%s]' % csi.curation.TF_instance.protein_accession) # TF 
        fields.append('%d' % csi.site_instance.start)                   # site start
        fields.append('%d' % csi.site_instance.end)                     # site end
        fields.append('%d' % csi.site_instance.strand)                  # site strand
        fields.append(csi.site_instance.seq)                     # site sequence
        regulated_genes = [reg.gene.locus_tag for reg in csi.regulation_set.all() if reg.evidence_type=="exp_verified"]
        fields.append(', '.join(regulated_genes) if regulated_genes else 'N/A')
        fields.append(csi.curation.publication.pmid)             # publication PMID
        fields.append(', '.join(tech.name for tech in csi.curation.experimental_techniques.all())) # tecniques
        fields.append('CollecTF:EXPSITE_%s' % utils.id2dbxref(int(csi.site_instance.site_id))),
        fields.append('collectf.umbc.edu' + reverse('browseapp.views.browse_by_site', kwargs={'dbxref_id':utils.id2dbxref(int(csi.site_instance.site_id))})),
        line = separator.join(fields)
        line = re.sub('[\n\r]', ' ', line) # remove newlines

        tsv_file.write(line.encode('utf-8'))
        tsv_file.write('\n')
        
    tsv_file.close()
        
    return tsv_file
