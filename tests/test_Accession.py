import unittest
from Accession import Accession
from pathlib import Path


class TestAccession(unittest.TestCase):

    accession_dict = {
        'NP_268346.1': 23, 'Q9CDN0.1': 33, 'pl_13.1': 33, 'sp|Q5XEC5|RL14_STRP6': 34, 'NP_1.1': 16, 'NP_6.2': 16,
        'NP_2.2': 25, 'NP_3.2': 42, 'sp|A0R8I7|RL16_BACAH': 42, 'sp|Q1GAM2|RNZ_LACDA': 43, 'NP_4.1': 27,
        'NP_5.2': 1911136, 'NP_7.1': 20, 'sp|Q97SU6|RL18_STRPN': 20, 'NP_8.2': 21, 'NP_9.2': 7137,
        'pl_12': 29,
    }

    def setUp(self):
        pass

    def test_read_accessions_gzipped(self):
        protaccession_path = Path.cwd() / 'data' / 'prot.accession2taxid_test.gz'
        pdbaccession_path = Path.cwd() / 'data' / 'pdb.accession2taxid_test.gz'
        accession = Accession({21, 42})
        accession.read_accessions(protaccession_path, pdbaccession_path, threads=4)
        self.assertEqual({'NP_3.2', 'sp|A0R8I7|RL16_BACAH', 'NP_8.2'}, accession.accessionIDs)

    def test_read_accessions_unzipped(self):
        protaccession_path = Path.cwd() / 'data' / 'prot.accession2taxid_test'
        pdbaccession_path = Path.cwd() / 'data' / 'pdb.accession2taxid_test.gz'
        accession = Accession({21, 42})
        accession.read_accessions(protaccession_path, pdbaccession_path, threads=4)
        self.assertEqual({'NP_3.2', 'sp|A0R8I7|RL16_BACAH', 'NP_8.2'}, accession.accessionIDs)

    def test_emptyDatabase(self):
        protaccession_path = Path.cwd() / 'data' / 'empty_accession.gz'
        pdbaccession_path = Path.cwd() / 'data' / 'empty_accession.gz'
        accession = Accession({21, 42})
        accession.read_accessions(protaccession_path, pdbaccession_path, threads=4)
        self.assertEqual(set(), accession.accessionIDs)
        protaccession_path = Path.cwd() / 'data' / 'empty_accession'
        accession = Accession({21, 42})
        accession.read_accessions(protaccession_path, pdbaccession_path, threads=4)
        self.assertEqual(set(), accession.accessionIDs)

        pass

if __name__ == '__main__':
    unittest.main()
