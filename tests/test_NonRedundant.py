import unittest
from NonRedundant import NonRedundant
from pathlib import Path
import shutil

class TestNonRedundant(unittest.TestCase):

    def setUp(self):
        pass

    def test_sequence_cleaner(self):
        test_db = Path.cwd() / 'data/redundant_test2.fasta'
        shutil.copyfile('data/redundant_test.fasta', 'data/redundant_test2.fasta')
        test_nr = 'data/redundant_test2_nr.fasta'
        file = open('data/non_redundant_test.fasta', "r")
        result_database = file.read()
        file.close()
        NonRedundant.sequence_cleaner(test_db, test_nr, False)
        file = open(test_nr)
        test_database = file.read()
        file.close()
        self.assertEqual(result_database.count('>'), test_database.count('>'))
        self.assertEqual(result_database.count('\n'), test_database.count('\n'))


if __name__ == '__main__':
    unittest.main()