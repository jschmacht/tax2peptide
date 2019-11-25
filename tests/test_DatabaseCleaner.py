import unittest
from DatabaseCleaner import DatabaseCleaner
from pathlib import Path
import shutil

class TestDatabaseCleaner(unittest.TestCase):

    def setUp(self):
        pass

    def test_non_redundant(self):
        test_db = Path.cwd() / 'data/redundant_test.fasta'
        file = open('data/non_redundant_test.fasta', "r")
        result_database = file.read()
        file.close()
        DatabaseCleaner.non_redundant(test_db, False)
        file = open(str(Path.cwd() / 'data/redundant_test_nr.fasta'))
        test_database = file.read()
        file.close()
        self.assertEqual(result_database.count('>'), test_database.count('>'))
        self.assertEqual(result_database.count('\n'), test_database.count('\n'))

    def test_reduce_ncbi_header(self):
        # copy example database
        DatabaseCleaner.reduce_header("./data/example_database_ncbi_rh2.fasta")
        result_database = open("./data/example_database_ncbi_rh.fasta", "r")
        test_database = open("./data/example_database_ncbi_reduced_header.fasta", "r")
        td = test_database.read()
        test_database.close()
        rd = result_database.read()
        result_database.close()
        self.assertEqual(td, rd)

if __name__ == '__main__':
    unittest.main()