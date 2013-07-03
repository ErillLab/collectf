# Create your views here.
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render
from django.shortcuts import render_to_response
from django.template import RequestContext
import models
from forms import ExportForm
import utils

@user_passes_test(lambda u: u.is_staff)
def export_tbl_view(request):
    if request.method=="GET":
        form = ExportForm()
        return render(request,
                      'ncbi_export.html',
                      {'form': form},
                      context_instance=RequestContext(request))
    else:
        form = ExportForm(request.POST)
        if form.is_valid():
            print form.cleaned_data
            export_tbl(form.cleaned_data['TF_instances'],
                       form.cleaned_data['genomes'])
            return render_to_response('ncbi_export_result.html',
                                      {'tbl': open('test.tbl').read()},
                                      context_instance=RequestContext(request))


def export_tbl(TF_instance, genome):
    """For given TF instance and genome, return all site instances in .tbl format"""

    # open tbl file to write
    tbl_file = open("test.tbl", 'w')
    tbl_file.write('>Feature test\n')
    
    # get all curation_site_instances
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome=genome,
        curation__TF_instance=TF_instance,
        curation__NCBI_submission_ready=True)
    # group curation_site_instance objects by site_instance
    site_instances = list(set(csi.site_instance for csi in curation_site_instances))
    for site_instance in site_instances:
        start, end = site_instance.start+1, site_instance.end+1
        if site_instance.strand == -1:
            start,end = end,start
        #tbl_file.write('%d %s\n' % (site_instance.strand, site_instance.seq))
        tbl_file.write('%d\t%d\tprotein_bind' % (start, end) + '\n')
        # all curation_site_instance objects of this site instance
        csis = [csi for csi in curation_site_instances if csi.site_instance==site_instance]
        # TF name
        assert all(csis[i].curation.TF.name == csis[0].curation.TF.name for i in xrange(len(csis)))
        tbl_file.write('\t\t\tbound_moiety\t%s\n' % (csis[0].curation.TF.name))
        tbl_file.write('\t\t\tnote\tTranscription factor binding site\n')
        # write experimental evidences
        for csi in csis:
            techs = csi.curation.experimental_techniques.all()
            tbl_file.write('\t\t\texperiment\t%s [PMID:%s]\n' % (', '.join(map(lambda t: t.name, techs)),
                                                                 csi.curation.publication.pmid))
        # write dbxref
        tbl_file.write('\t\t\tdb_xref\tCollecTF:EXPSITE_%s\n' % utils.id2dbxref(int(site_instance.site_id)))
        # write regulation note
        evidence4regulation = set([reg.gene.locus_tag for csi in csis for reg in csi.regulation_set.all() if reg.evidence_type=="exp_verified"])
        if evidence4regulation:
            tbl_file.write('\t\t\tnote\tEvidence of regulation for: %s\n' % (', '.join(evidence4regulation)))
                       

    tbl_file.close()
