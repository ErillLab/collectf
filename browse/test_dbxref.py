from django.test import TestCase

from browse import dbxref

class UniprotDbxrefTest(TestCase):

    def test_converts_id_to_uniprot_dbxref(self):
        self.assertEqual(dbxref.to_uniprot_dbxref(1), 'EXPREG_00000010')
        self.assertEqual(dbxref.to_uniprot_dbxref(100), 'EXPREG_00000640')

    def test_converts_uniprot_dbxref_to_id(self):
        self.assertEqual(dbxref.from_uniprot_dbxref('00000010'), 1)
        self.assertEqual(dbxref.from_uniprot_dbxref('00000640'), 100)
