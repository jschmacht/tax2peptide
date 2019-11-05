import unittest
from WriteCustomDB import WriteCustomDB


class TestCustomTaxonDB(unittest.TestCase):

    def setUp(self):
        pass

    def test_emptyTaxIDs(self):
        taxIDs=set()
        write_db = WriteCustomDB("./data/example_database_uniprot.fasta", "./data/customized_uniprot.fasta", taxIDs)
        write_db.read_database(True, gzipped=False, threads=3)
        file = open("./data/customized_uniprot.fasta", "r")
        self.assertEqual(file.read(), "")
        file.close()

    def test_emptyDatabase(self):
        taxIDs = {5, 6, 26}
        write_db = WriteCustomDB("./data/empty_database.fasta", "./data/customized_uniprot.fasta", taxIDs)
        write_db.read_database(True, gzipped=False, threads=3)
        file = open("./data/customized_uniprot.fasta", "r")
        self.assertEqual(file.read(), "")
        file.close()

    def test_ncbi_db_search(self):
        accessionIDs = {'NP_268346.1', 'NP_5.2', 'NP_3.2', 'AAN62213.1', 'NP_9.2', 'pl_12'}
        #accessionIDs = {'AAN62213.1'}
        write_db = WriteCustomDB("./data/example_database_ncbi.fasta", "./data/customized_ncbi.fasta")
        write_db.read_database(with_taxon_ID=False, gzipped=False, threads=3, accessions=accessionIDs)
        test_database = open("./data/result_test_ncbi.fasta", "r")
        td = test_database.read()
        test_database.close()
        written_database = open("./data/customized_ncbi.fasta", "r")
        wd = written_database.read()
        written_database.close()
        self.assertEqual(wd.count('>'), td.count('>'))
        self.assertEqual(wd.count('\n'), td.count('\n'))

    def test_ncbi_db_search_gzipped(self):
        accessionIDs = {'NP_268346.1', 'NP_5.2', 'NP_3.2', 'AAN62213.1', 'NP_9.2', 'pl_12'}
        write_db = WriteCustomDB("./data/example_database_ncbi.fasta.gz", "./data/customized_ncbi.fasta")
        write_db.read_database(False, threads=2, accessions=accessionIDs, gzipped=True)
        test_database = open("./data/result_test_ncbi.fasta", "r")
        td = test_database.read()
        written_database = open("./data/customized_ncbi.fasta", "r")
        wd = written_database.read()
        self.assertEqual(wd.count('>'), td.count('>'))
        self.assertEqual(wd.count('\n'), td.count('\n'))
        test_database.close()
        written_database.close()

    def test_uniprot_db_search(self):
        taxIDs = {6, 7, 34, 26}
        write_db = WriteCustomDB("./data/example_database_uniprot.fasta", "./data/customized_uniprot.fasta", taxIDs)
        write_db.read_database(True, gzipped=False)
        test_database = open("./data/result_test_uniprot.fasta", "r")
        td = test_database.read()
        test_database.close()
        written_database = open("./data/customized_uniprot.fasta", "r")
        wd = written_database.read()
        written_database.close()
        self.assertEqual(wd.count('>'), td.count('>'))
        self.assertEqual(wd.count('\n'), td.count('\n'))

    def test_uniprot_db_search_gzipped(self):
        taxIDs = {6, 7, 34, 26}
        write_db = WriteCustomDB("./data/example_database_uniprot.fasta.gz", "./data/customized_uniprot.fasta", taxIDs)
        write_db.read_database('uniprot', gzipped=True, threads=3)
        test_database = open("./data/result_test_uniprot.fasta", "r")
        td = test_database.read()
        written_database = open("./data/customized_uniprot.fasta", "r")
        wd = written_database.read()
        self.assertEqual(wd.count('>'), td.count('>'))
        self.assertEqual(wd.count('\n'), td.count('\n'))
        test_database.close()
        written_database.close()


if __name__ == '__main__':
    unittest.main()
