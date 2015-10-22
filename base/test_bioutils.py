from django.test import TestCase
from base import bioutils

class SequenceTest(TestCase):

    def test_performs_reverse_complement(self):
        seq = 'ACGTAAA'
        self.assertEqual(bioutils.reverse_complement(seq), 'TTTACGT')

    def test_reverse_complements_empty_sequence(self):
        self.assertEqual(bioutils.reverse_complement(''), '')


class EntrezTest(TestCase):
    def test_fetches_pubmed_record(self):
        record = bioutils.get_pubmed(1)
        self.assertEqual(
            record['Title'],
            'Formate assay in body fluids: application in methanol poisoning.')

    def test_none_on_invalid_pmid(self):
        record = bioutils.get_pubmed(-1)
        self.assertIsNone(record)
    
