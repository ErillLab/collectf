from collectfapp import models
from collectfapp import bioutils
from collectfapp import sitesearch


def site_fix():
    same = 0
    not_same = 0
    for csi in models.Curation_SiteInstance.objects.all().iterator():
        #print csi.site_instance.site_id
        seq_in_genome = csi.site_instance.genome.sequence[csi.site_instance.start:csi.site_instance.end+1]
        if csi.site_instance.strand == -1:
            seq_in_genome = bioutils.reverse_complement(seq_in_genome)

        if not seq_in_genome == csi.site_instance.seq:
            #print csi.site_instance.seq, csi.site_instance.start, csi.site_instance.end, seq_in_genome
            print csi.site_instance.site_id
            start = csi.site_instance.start-1
            end= csi.site_instance.end-2
            new_seq_in_genome = csi.site_instance.genome.sequence[start:end+1]
            if csi.site_instance.strand==-1:
                new_seq_in_genome = bioutils.reverse_complement(new_seq_in_genome)
            assert new_seq_in_genome == csi.site_instance.seq, new_seq_in_genome + ' | ' + csi.site_instance.seq
            print 'fixing', csi.site_instance.site_id
            csi.site_instance.start = start
            csi.site_instance.end = end
            csi.site_instance.save()
            

        """
        if csi.site_instance.seq != seq_in_genome:
            print csi.site_instance.genome.strain.name
            print csi.site_instance.seq
            print seq_in_genome, csi.site_instance.strand
            print csi.annotated_seq
            print ''
            
            #print 'change?'
            #raw_input()
            #csi.site_instance.end -= 1
            #csi.site_instance.save()
            #csi.annotated_seq, csi.site_instance.seq = csi.site_instance.seq, csi.annotated_seq
            #csi.site_instance.save()
            #csi.save()
        
        if csi.site_instance.seq == csi.annotated_seq:
            same += 1
        else:
            not_same += 1
    print 'same', same, 'not_same', not_same
    """


def check_regulation_duplicates():
    i = 0
    for reg in models.Regulation.objects.all():
        if i % 100 == 0: print i
        i += 1
        same_genes = models.Gene.objects.filter(gene_accession=reg.gene.gene_accession)

        if len(same_genes) > 1:
            min_gene_id = min(g.gene_id for g in same_genes)
            print min_gene_id, reg.gene.gene_id

def check_multiple_genes():
    unique_gene_accessions = list(set(g.gene_accession for g in models.Gene.objects.all()))
    print len(models.Gene.objects.all().iterator())
    print len(unique_gene_accessions)

    d = dict()
    for g in models.Gene.objects.all().iterator():
        d[g.gene_accession] = d.get(g.gene_accession, []) + [g]

    print 'doing'
    
    for i,ga in enumerate(unique_gene_accessions):

        #print 'processing %d/%d' % (i, len(unique_gene_accessions))
        all_copies = d[ga]
        if len(all_copies) > 1:


            min_gene_id = min(g.gene_id for g in all_copies)

            to_delete = [x for x in all_copies if x.gene_id != min_gene_id]
            
#            for x in to_delete:
#                print x
#                x.delete()

            
def check_multiple_genomes():
    unique_genome_accessions = list(set(g.genome_accession.split('.')[0] for g in models.Genome.objects.all()))
    print len(models.Genome.objects.all())
    print len(unique_genome_accessions)
        
def run():
    #check_multiple_genes()
    #check_multiple_genomes()
    site_fix()
    print "I am a script"
    
