from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import time
import os
import pickle
from collectf import settings
from base64 import b64encode
import lasagna



Entrez.email = "sefakilic@gmail.com"

def get_pubmed(pmid):
    """Retrieve pubmed publication from NCBI database."""
    try:
        handle = Entrez.esummary(db="pubmed", id=pmid)
        record = Entrez.read(handle)
        return record[0]
    except RuntimeError:
        print 'err'
        return None

def get_genome(accession):
    """Retrieve genome from NCBI database"""
    try:
        print accession
        h = Entrez.efetch(db='nuccore', id=accession, retmode='gbwithparts', rettype='text')
        seq_record = SeqIO.read(h, 'gb')
        h.close()
        return seq_record
    except:
        print 'Something went wrong during refseq retrieval'
        return None

def get_TF(accession):
    """Retrieve transcription factor from NCBI database"""
    try:
        h = Entrez.efetch(db='protein', id=accession, retmode='gb', rettype='text')
        seq_record = SeqIO.read(h, 'gb')
        h.close()
        return seq_record
    except:
        print 'Something went wrong during TF sequence retrieval'
        return None


def get_gene_id(feature):
    # extract gene id
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
    # Occasionally, when collecTF tries to retrieve all gene summaries, NCBI refuses
    # to return all them. There should be some sort of limit for a query. Therefore,
    # the gene list summary query is chunked into 1000 gene pieces.
    
    # edit: After several different ways to fix it, the best way is breaking up the
    # list into batches of size 1000 AND retry if any batch fails for some reason.
    runtime_error = 0  # number of runtime errors during Entrez esummary
    
    while runtime_error < 10:  # if runtime error is consistent, there is no point trying again
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

    assert runtime_error == 0
    return records

def get_genes(genome_rec):
    """Get list of all genes"""
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
    """Given genome record from NCBI db, get organism name"""
    return genome_record.annotations['organism']

def get_org_taxon_depr(genome_record):
    """Given genome rec from NCBI db, get organism taxonomy id"""
    org = get_org_name(genome_record)
    org = org.replace('(', ' ')
    org = org.replace(')', ' ')
    handle = Entrez.esearch(db='taxonomy', term=org)
    rec = Entrez.read(handle)
    assert int(rec['Count']) == 1
    return rec['IdList'][0]

def get_org_taxon(genome_record):
    """Given genome record, find organism taxonomy id using Elink utility"""
    try:
        gi = genome_record.annotations['gi']
        r = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='nuccore', id=gi, linkname='nuccore_taxonomy'))
        assert len(r) == 1
        tax_id = r[0]['LinkSetDb'][0]['Link'][0]['Id']
    except:
        tax_id = None
    return tax_id

def get_taxon_info(tax_id):
    try:
        handler = Entrez.efetch(db="Taxonomy", id=str(tax_id), retmode="xml")
        records = Entrez.read(handler)
        record = records[0]
        lineage = record['LineageEx']
        order = [x for x in lineage if x['Rank']=='order']
        order_name = order[0]['ScientificName']
    except:
        order_name = 'other'
    print tax_id, order_name
    return order_name


TAXON_FILE = os.path.join(settings.PICKLE_ROOT, "taxon.pickle")

def get_all_taxon_info(tax_ids):
    """For all taxonomy ids in tax_ids, get taxonomy id information and pickle it"""
    import pickle
    taxon = {}
    try:
        taxon = pickle.load(open(TAXON_FILE))
    except:
        pass
    for tax_id in [tax_id for tax_id in tax_ids if tax_id not in taxon]:
        taxon[tax_id] = get_taxon_info(tax_id)
    pickle.dump(taxon, open(TAXON_FILE, 'w'))
    return taxon

def get_taxon_info_from_file(tax_id):
    try:
        taxon = pickle.load(open(TAXON_FILE))
        if tax_id not in taxon:
            raise Exception("tax id not found")
    except:
        taxon = get_all_taxon_info([tax_id])
    return taxon[tax_id]
    
def to_fasta(seqs):
    """
    FASTA representation of motif
    """
    str = ""
    for i,inst in enumerate(seqs):
        str = str + ">instance%d\n"%i + inst + "\n"
    return str   

def weblogo(sequences, format="PNG"):
    """
    uses the Berkeley weblogo service to download and save a weblogo of itself
    
    requires an internet connection.
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
    
    #for k,v in kwds.iteritems():
    #    values[k]=str(v)

    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    im=response.read()

    return im
    
def weblogo_uri(sequences):
    """Generate the weblogo and make it ready for direct embed into response html"""
    image_data = weblogo(sequences)
    encoded = b64encode(image_data)
    mime = "image/png"
    
    return ("data:" + mime + ';' + "base64," + encoded)


def get_overlap(loca, locb):
    """Given two regions, return the length of overlap of them"""
    return max(0, min(loca[1], locb[1]) - max(loca[0], locb[0]))

def overlap_site_meta_site(curation_site_instance, meta_site_instance, overlap_th=0.8):
    """Given a site instance (site_instance) and a list of site instances
    (meta_site_instance), return whether site instance overlaps enough with any site
    instance in the meta_site_instance."""

    def location(curation_site_instance):
        return (curation_site_instance.site_instance.start,
                curation_site_instance.site_instance.end)
    
    return any(get_overlap(location(curation_site_instance), location(ms)) for ms in meta_site_instance)
