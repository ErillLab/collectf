from django.test import TestCase
from base import bioutils

class SequenceTest(TestCase):
    def test_performs_reverse_complement(self):
        seq = 'ACGTAAA'
        self.assertEqual(bioutils.reverse_complement(seq), 'TTTACGT')

    def test_reverse_complements_empty_sequence(self):
        self.assertEqual(bioutils.reverse_complement(''), '')

    def test_converts_to_FASTA_format(self):
        seqs = ['ACGT', 'ACCCA']
        self.assertEqual(bioutils.to_fasta(seqs),
                         '>instance0\nACGT\n>instance1\nACCCA\n')

class EntrezTest(TestCase):
    def setUp(self):
        self.genome_record = bioutils.get_genome('NC_000913.3')
        
    def test_fetches_pubmed_record(self):
        record = bioutils.get_pubmed(1)
        self.assertEqual(
            record['Title'],
            'Formate assay in body fluids: application in methanol poisoning.')

    def test_returns_nothing_on_invalid_pmid(self):
        record = bioutils.get_pubmed(-1)
        self.assertIsNone(record)
    
    def test_fetches_genome_record(self):
        self.assertEqual(self.genome_record.id, 'NC_000913.3')
        self.assertEqual(str(self.genome_record.seq[:5]), 'AGCTT')

    def test_returns_nothing_on_invalid_genome_accession(self):
        record = bioutils.get_genome('NC_invalid')
        self.assertIsNone(record)

    def test_fetches_TF_record(self):
        accession = 'NP_233336.1'
        record = bioutils.get_TF(accession)
        self.assertEqual(record.id, accession)
        self.assertEqual(record.seq[:5], 'MKDEN')

    def test_returns_none_on_invalid_TF_accession(self):
        record = bioutils.get_TF('invalid_TF_accession')
        self.assertIsNone(record)

    def test_fetches_TF_uniprot_record(self):
        record = bioutils.get_uniprot_TF('Q9KKZ8')
        self.assertIsInstance(record, unicode)
        self.assertIn('LuxR family', record)

    def test_returns_empty_on_invalid_uniprot_accession(self):
        record = bioutils.get_uniprot_TF('invalid_uniprot_acc')
        self.assertIsInstance(record, unicode)
        self.assertEqual(record, '')

    def test_gets_gene_id_from_feature(self):
        gene_features = [feature for feature in self.genome_record.features
                         if feature.type == 'gene']
        gene_id = bioutils.get_gene_id(gene_features[0])
        self.assertEqual(gene_id, '944742')

    def test_gets_genes(self):
        genes = bioutils.get_genes(self.genome_record)
        self.assertEqual(len(genes), 4498)
        expected = {'locus_tag': 'b0001',
                    'name': 'thrL',
                    'gene_accession': '944742',
                    'description': 'thr operon leader peptide',
                    'start': 189,
                    'end': 255,
                    'strand': 1}
        self.assertEqual(genes[0], expected)

    def test_gets_org_name(self):
        self.assertEqual(bioutils.get_org_name(self.genome_record),
                         'Escherichia coli str. K-12 substr. MG1655')

    def test_gets_org_taxon(self):
        self.assertEqual(bioutils.get_organism_taxon(self.genome_record),
                         '511145')

    def test_gets_org_taxon_from_TF_accession(self):
        acc = 'WP_010970539.1'
        self.assertEqual(bioutils.TF_accession_to_org_taxon(acc), '382')

    
