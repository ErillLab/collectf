# Create your views here.
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render
from django.shortcuts import render_to_response
from django.template import RequestContext
import models
import os
from forms import ExportForm
from baseapp import utils
import tempfile

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
            tbl_data = export_tbl(form.cleaned_data['TF_instances'], form.cleaned_data['genomes'])
            return render_to_response('ncbi_export_result.html',
                                      {'tbl': tbl_data,
                                       'TF_instance': form.cleaned_data['TF_instances'],
                                       'genome': form.cleaned_data['genomes'],
                                       },
                                      context_instance=RequestContext(request))


def export_tbl(TF_instance, genome):
    """For given TF instance and genome, return all site instances in .tbl format"""

    # open tbl file to write
    tbl_file = tempfile.TemporaryFile()
    tbl_file.write('>Feature prot_%s_genome_%s\n' % (TF_instance.protein_accession, genome.genome_accession))
    
    # get all curation_site_instances
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome=genome,
        curation__TF_instance=TF_instance,
        curation__NCBI_submission_ready=True,
        is_motif_associated=True,
        curation__experimental_techniques__preset_function__in=['binding', 'expression'])

    # group curation_site_instance objects by meta-site-instance
    meta_site_instances = models.MetaSiteInstance.objects.filter(curation_siteinstance__in=curation_site_instances).distinct().order_by('start')
    for meta_site_instance in meta_site_instances:
        # all curation_site_instance objects of this meta_site_instance
        csis = meta_site_instance.curation_siteinstance_set.all()
        assert csis.count() > 0
        start = min(map(lambda csi: csi.site_instance.start, csis))+1
        end = max(map(lambda csi: csi.site_instance.end, csis))+1
        # TF name consistency check
        if any (csis[i].curation.TF != csis[0].curation.TF for i in xrange(len(csis))):
            tbl_file.truncate() # remove the contents (if any)
            tbl_file.write("Inconsistent TF - TF_instance matches. This TF_instance is linked to more than one TF.\n")
            return tbl_file.read()
        tbl_file.write('%d\t%d\tprotein_bind' % (start, end) + '\n')
        tbl_file.write('\t\t\tbound_moiety\t%s\n' % (csis[0].curation.TF.name))
        tbl_file.write('\t\t\tnote\tTranscription factor binding site\n')
        # write experimental evidences
        experiments = {}
        for exp in models.ExperimentalTechnique.objects.filter(preset_function__in=['binding', 'expression']):
            filtered_csis = [csi for csi in csis if exp in csi.curation.experimental_techniques.all()]
            experiments[exp] = list(set([csi.curation.publication.pmid for csi in filtered_csis]))
        for exp,pmids in experiments.items():
            if not pmids: continue
            tbl_file.write('\t\t\texperiment\t%s [PMID: %s]\n' % (exp.name, ', '.join(pmids)))

        # write regulation note
        evidence4regulation = set([reg.gene.locus_tag for csi in csis for reg in csi.regulation_set.all() if reg.evidence_type=="exp_verified"])
        if evidence4regulation:
            tbl_file.write('\t\t\tnote\tEvidence of regulation for: %s\n' % (', '.join(evidence4regulation)))

        # write dbxref
        tbl_file.write('\t\t\tdb_xref\t%s\n' % utils.id2dbxref(int(meta_site_instance.pk)))

    tbl_file.seek(0) # goto beginnning of the file
    return tbl_file.read()

