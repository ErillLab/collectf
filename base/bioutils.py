"""This module contains generic functions that are used in collectf apps, for
processing sequence information. Most of them are built on top of Biopython."""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from base64 import b64encode
import models
import lasagna
from subprocess import Popen, PIPE, STDOUT

Entrez.email = "sefa1@umbc.edu"

def reverse_complement(seq):
    """Reverse complement of a sequence"""
    return Seq(seq).reverse_complement().tostring()

def get_pubmed(pmid):
    """Retrieve pubmed publication from NCBI database."""
    try:
        handle = Entrez.esummary(db="pubmed", id=pmid)
        record = Entrez.read(handle)
        return record[0]
    except RuntimeError:
        return None

def get_genome(accession):
    """Retrieve genome record from NCBI database."""
    try:
        h = Entrez.efetch(db='nuccore', id=accession, retmode='gbwithparts', rettype='text')
        seq_record = SeqIO.read(h, 'gb')
        h.close()
        return seq_record
    except:
        return None

def get_TF(accession):
    """Retrieve transcription factor from NCBI database."""
    try:
        h = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='gb')
        seq_record = SeqIO.read(h, 'gb')
        h.close()
        return seq_record
    except RuntimeError as e:
        print e
        return None

def get_gene_id(feature):
    """Given a Biopython SeqFeature object, extract the gene id."""
    i = 0
    while not feature.qualifiers['db_xref'][i].startswith('GeneID:'):
        i += 1
    gene_id = feature.qualifiers['db_xref'][i][7:]
    return gene_id

def get_gene_annotation(id_list):
    """Annotates Entrez Gene ids using Bio.Entrez, in particular epost (to submit the
    data to NCBI) and esummary to retrieve the information. Returns a list gene
    summary objects."""
    epost_result = Entrez.read(Entrez.epost("gene", id=','.join(id_list)))

    # Occasionally, when CollecTF tries to retrieve all gene summaries, NCBI
    # refuses to return all them. There must be some sort of limit for a
    # query. Therefore, the gene list summary query is chunked into 1000 gene
    # pieces.
    # Edit: After several different ways to fix it, the best way is breaking up
    # the list into batches of size 1000 AND retry if any batch fails for some
    # reason.
    runtime_error = 0  # number of runtime errors during Entrez esummary
    while runtime_error < 10:  # if runtime error is consistent, there is no
                               # point trying again
        try:
            request = Entrez.esummary(db="gene", webenv=epost_result["WebEnv"],
                                      query_key=epost_result["QueryKey"])
            records = Entrez.read(request)
        except RuntimeError as e:
            print "Error occurred during epost+esummary:", e
            print "Trying again."
            runtime_error += 1
        else:
            runtime_error = 0
            break

    assert runtime_error == 0  # make sure it is completed successfully
    return records

def get_genes(genome_rec):
    """Given a genome record object, get list of all genes."""
    genes = [] # return list of genes
    # Use Entrez post method, because get method has limitation on url length
    # use Epost to post list of ids first
    # get gene ids
    gene_features = [f for f in genome_rec.features if f.type == 'gene']
    gids = [get_gene_id(f) for f in gene_features]
    
    recs = []
    chunk_size = 1000
    for start in xrange(0, len(gids), chunk_size):
        end = min(len(gids), start+chunk_size)
        #print gids[start:end]
        recs = recs + get_gene_annotation(gids[start:end])
        
    # two sources of gene data: Entrez epost and gene features from genome record
    for gid, feat, rec in zip(gids, gene_features, recs):
        assert rec['Id'] == gid
        genes.append({'gene_accession': gid,
                      'name': rec['Name'],
                      'description': rec['Description'],
                      'start': feat.location.start.position,
                      'end': feat.location.end.position,
                      'strand': feat.strand,
                      'locus_tag': ','.join(feat.qualifiers['locus_tag'])})
    return genes

def get_org_name(genome_record):
    """Given genome record from NCBI db, get organism name."""
    return genome_record.annotations['organism']

def get_organism_taxon(genome_record):
    """Given genome record, find organism taxonomy id using Elink utility."""
    try:
        # Get taxonomy id using Entrez elink utility
        gi = genome_record.annotations['gi']
        r = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='nuccore',
                                     id=gi, linkname='nuccore_taxonomy'))
        assert len(r) == 1
        tax_id = r[0]['LinkSetDb'][0]['Link'][0]['Id']
    except:
        tax_id = None
    return tax_id

def TF_accession_to_org_taxon(TF_accession):
    """Given a TF accession number, find the organism taxonomy id using Elink utility"""
    TF_record = get_TF(TF_accession)
    if not TF_record:
        raise RuntimeError # raise exception if the TF record can not be fetched
    # Get taxonomy id using Entrez elink utility
    gi = TF_record.annotations['gi']
    r = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='protein',
                    id=gi, linkname='protein_taxonomy'))
    assert len(r) == 1
    tax_id = r[0]['LinkSetDb'][0]['Link'][0]['Id']
    return tax_id

def to_fasta(seqs):
    """Given a collection of sequences, put them into the FASTA format. It is
    used for weblogo generation."""
    str = ""
    for i,inst in enumerate(seqs):
        str = str + ">instance%d\n"%i + inst + "\n"
    return str

def weblogo(sequences):
    """Given a collection of site sequences, generate the sequence logo, using
    weblogo program that is locally installed."""
    al = to_fasta(sequences)    
    p = Popen(['/usr/bin/weblogo', '-F', 'png', '-s', 'LARGE', '-c', 'classic', '--errorbars', 'YES'],
              stdout=PIPE, stdin=PIPE, stderr=PIPE, close_fds=True)
    stdout_data, stderr_data = p.communicate(input=al)
    return stdout_data

def weblogo_deprecated(sequences, format="PNG"):
    """
    Weblogo generation (deprecated).
    
    It uses the Berkeley weblogo service to download and save a weblogo of
    itself, requires an internet connection. It has been replaced with the
    system call that runs weblogo locally.
    
    The parameters from **kwds are passed directly to the weblogo server.
    """
    import urllib
    import urllib2

    # if all sequences don't have the same length, apply LASAGNA
    # use LASAGNA to align sites
    aligned, idxAligned, strands = lasagna.LASAGNA(map(lambda s: str(s.lower()), sequences), 0)
    trimmed = lasagna.TrimAlignment(aligned) if len(aligned) > 1 else aligned
    trimmed = [s.upper() for s in trimmed]
    
    assert all(len(seq) == len(trimmed[0]) for seq in trimmed), "sequences do not have the same length"

    al = to_fasta(trimmed)
    
    url = 'http://weblogo.berkeley.edu/logo.cgi'
    values = {'sequence' : al,
              'format' : format,
              'logowidth' : '18',
              'logoheight' : '5',
              'logounits' : 'cm',
              'kind' : 'AUTO',
              'firstnum' : "1",
              'command' : 'Create Logo',
              'smallsamplecorrection' : "on",
              'symbolsperline' : 32,
              'res' : '96',
              'res_units' : 'ppi',
              'antialias' : 'on',
              'title' : '',
              'barbits' : '',
              'xaxis': 'on',
              'xaxis_label'  : '',
              'yaxis': 'on',
              'yaxis_label' : '',
              'showends' : 'on',
              'shrink' : '0.5',
              'fineprint' : 'on',
              'ticbits' : '1',
              'colorscheme' : 'DEFAULT',
              'color1' : 'green',
              'color2' : 'blue',
              'color3' : 'red',
              'color4' : 'black',
              'color5' : 'purple',
              'color6' : 'orange',
              'color1' : 'black',
              }
    
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    im=response.read()
    return im

def weblogo_uri(sequences):
    """Generate the weblogo and make it ready for direct embed into response
    HTML."""
    image_data = weblogo_deprecated(sequences)
    encoded = b64encode(image_data)
    mime = "image/png"
    return ("data:" + mime + ';' + "base64," + encoded)

def run_lasagna(site_instances, trim=True):
    """Given a list of models.SiteInstance objects, run LASAGNA algorithm and return the
    aligned sequences. Here is the link to the paper:
    http://www.biomedcentral.com/1471-2105/14/108
    """


    sites = [x.seq_lower for x in site_instances]
    aligned, idxAligned, strands = lasagna.LASAGNA(sites, 0)
    aligned = [s.upper() for s in aligned]

    # if you asked for trimmed alignment (only the region that has no gaps)
    if not trim: return aligned

    # If asked for non-trimmed, the gaps in the alignment should be recovered
    # using genome sequence
    assert map(int, idxAligned) == range(len(idxAligned)), "lasagna sites are not sorted"
    # batch fetch of genomes from the database
    queryset = models.SiteInstance.objects.filter(pk__in=[x.pk for x in site_instances])
    site_genome_dict = dict(v for v in queryset.values_list("site_id", "genome").distinct())
    recovered = [fill_gaps(site_instance, aligned_site, aligned_strand)
                 for (site_instance, aligned_site, aligned_strand) in
                 zip(site_instances, aligned, strands)]

    return recovered



def extend_site(site_instance, genome_seq, n=250):
    """Extend the site instance sequence by n bases both sides."""
    seq = genome_seq[site_instance.start-n: site_instance.end+n+1]
    return (seq if site_instance.strand==1 else reverse_complement(seq))

def fill_gaps(site_instance, aligned_site, aligned_strand):
    """Fill the gaps in the alignment."""
    genome_seq = site_instance.get_genome_sequence()
    seq = site_instance.seq
    original_seq = seq if aligned_strand=='+' else reverse_complement(seq)
    num_left_gaps = 0
    num_right_gaps = 0
    # count left and right gaps
    while aligned_site[num_left_gaps] == '-':
        num_left_gaps += 1
    while aligned_site[-num_right_gaps-1] == '-':
        num_right_gaps += 1
    # check if correct
    assert '-' * num_left_gaps + str(original_seq) + '-'*num_right_gaps == aligned_site
    # extend site
    n = max(num_left_gaps, num_right_gaps)
    # extend to both sides
    extended_site = extend_site(site_instance, genome_seq, n)
    if aligned_strand == '-':
        extended_site = reverse_complement(extended_site)

    if num_right_gaps == n:
        recovered_site = extended_site[n-num_left_gaps:]
    else:
        recovered_site = extended_site[n-num_left_gaps: num_right_gaps-n]
    # some checks
    assert len(recovered_site) == len(aligned_site), '%s %d' % (str(site_instance), n)
    assert recovered_site in extended_site
    return str(recovered_site)


def build_motif(seqs):
    """Create motif from sequences"""
    m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
    for seq in seqs:
        try:
            m.add_instance(Seq(seq, m.alphabet))
        except:
            print "Diff motif size length?"
            return None
    m.make_counts_from_instances()
    return m

def score_sequence(motif, sequence):
    """Given Biopython motif object and a sequence match, return PWM score of the match"""
    return motif.score_hit(sequence, position=0)


def degenerate_consensus(motif):
    """Grabbed and modified from Biopython modules library.  Following the rules
    adapted from D. R. Cavener: "Comparison of the consensus sequence flanking
    translational start sites in Drosophila and vertebrates."  Nucleic Acids
    Research 15(4): 1353-1361. (1987).  The same rules are used by TRANSFAC.
    """
    degenerate_nucleotide = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'AC': 'M',
        'AG': 'R',
        'AT': 'W',
        'CG': 'S',
        'CT': 'Y',
        'GT': 'K',
        'ACG': 'V',
        'ACT': 'H',
        'AGT': 'D',
        'CGT': 'B',
        'ACGT': 'N',
    }
    sequence = ""
    for i in range(motif.length):
        def get(nucleotide):
            return motif.counts[nucleotide][i]
        nucleotides = sorted(motif.counts, key=get, reverse=True)
        counts = [motif.counts[c][i] for c in nucleotides]
        # Follow the Cavener rules:
        if counts[0] >= sum(counts[1:]) and counts[0] >= 2*counts[1]:
            key = nucleotides[0]
        elif 4*sum(counts[:2]) > 3*sum(counts):
            key = "".join(sorted(nucleotides[:2]))
        elif counts[3]==0:
            key = "".join(sorted(nucleotides[:3]))
        else:
            key = "ACGT"
        nucleotide = degenerate_nucleotide[key]
        sequence += nucleotide
    return sequence
