import unittest
from TaxonGraph import TaxonGraph


class TestCreateGraph(unittest.TestCase):

    def setUp(self):
        pass

    def test_emptyFile(self):
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/emptyFile.tar.gz')
        self.assertEqual(taxon_graph.parent_child_graph, {})
        self.assertEqual(taxon_graph.child_rank_dict, {})
        self.assertEqual(taxon_graph.oldtax_newtax_dict, {})
        self.assertEqual(taxon_graph.taxon_name_dict, {})
        self.assertEqual(taxon_graph.child_parent_dict, {})

    def test_create_graph(self):
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        expected_child_parent_dict = {
            1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 3, 8: 3, 9: 3, 10: 4, 11: 4, 12: 4, 13: 4, 14: 5, 15: 6, 16: 15,
            17: 7, 18: 7, 19: 8, 20: 9, 21: 20, 22: 21, 23: 10, 24: 10, 25: 10, 26: 11, 27: 11, 28: 12, 29: 28, 30: 29,
            31: 13, 32: 13, 33: 13, 34: 33, 35: 34, 36: 34, 37: 2, 38: 14
        }
        expected_child_rank_dict = {
            1: 'no rank', 2: 'superkingdom', 3: 'genus', 4: 'genus', 5: 'family', 6: 'no rank', 7: 'section', 8: 'section',
            9: 'no rank', 10: 'species', 11: 'section', 12: 'species', 13: 'no rank', 14: 'species', 15: 'species',
            16: 'subspecies', 17: 'species', 18: 'species', 19: 'species', 20: 'no rank', 21: 'species',
            22: 'subspecies', 23: 'subspecies', 24: 'subspecies', 25: 'no rank', 26: 'species', 27: 'species',
            28: 'no rank', 29: 'subspecies', 30: 'no rank', 31: 'species', 32: 'subspecies', 33: 'species',
            34: 'subspecies', 35:'no rank', 36: 'no rank', 37: 'no rank', 38: 'subspecies'
        }
        expected_oldtax_newtax_dict = {45: 11, 36: 12, 100: 2, 44: 34}
        expected_parent_child_graph = {
            1: [2, 3, 4], 2: [5, 6, 37], 3: [7, 8, 9], 4: [10, 11, 12, 13], 5: [14], 6: [15], 7: [17, 18], 8: [19], 9: [20],
            10: [23, 24, 25], 11: [26, 27], 12: [28], 13: [31, 32, 33], 14: [38], 15: [16], 16: [], 17: [], 18: [], 19: [],
            20: [21], 21: [22], 22: [], 23: [], 24: [], 25: [], 26: [], 27: [], 28: [29], 29: [30], 30: [], 31: [],
            32: [], 33: [34], 34: [35, 36], 35: [], 36: [], 37: [], 38: []
        }
        expected_names_dict = {
            1: 'root', 2: 'Bacteria', 3: 'Stigmatella', 4: 'Archangiaceae', 5: 'Archangium disciforme', 6: 'Azorhizobium',
            7: 'Azorhizobium caulinodans', 8: 'Myxococcus macrosporus', 9: 'Buchnera aphidicola', 10: 'Cellvibrio',
            11: 'Cellulomonas gilvus', 12: 'Stigmatella', 13: 'Dictyoglomus', 14: 'Dictyoglomus thermophilum',
            15: 'Melittangium lichenicola', 16: 'Methylophilus', 17: 'Methylophilus methylotrophus', 18: 'Pelobacter',
            19: 'Pelobacter carbinolicus', 20: 'Phenylobacterium', 21:'Phenylobacterium immobile', 22: 'Shewanella',
            23: 'Shewanella colwelliana', 24: 'Shewanella putrefaciens', 25: 'Shewanella hanedai',
            26:'Stigmatella aurantiaca', 27: 'halophilic eubacterium NRCC 26227', 28: 'halophilic eubacterium',
            29: 'Myxococcales', 30: 'Cystobacter', 31: 'Myxococcaceae', 32: 'Myxococcus', 33: 'Myxococcus fulvus',
            34: 'Myxococcus xanthus', 35: 'Cystobacter fuscus', 36: 'Melittangium', 37: '37', 38: '38'
        }
        self.assertEqual(expected_parent_child_graph, taxon_graph.parent_child_graph)
        self.assertEqual(expected_child_rank_dict, taxon_graph.child_rank_dict)
        self.assertEqual(expected_oldtax_newtax_dict, taxon_graph.oldtax_newtax_dict)
        self.assertEqual(expected_names_dict, taxon_graph.taxon_name_dict)
        self.assertEqual(expected_child_parent_dict, taxon_graph.child_parent_dict)

    def test_find_species_taxIDs(self):
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        self.assertEqual({9, 20, 21, 22}, taxon_graph.find_taxIDs(9))

    def test_find_species_taxIDs_error(self):
        taxon_graph = TaxonGraph()
        with self.assertRaises(Exception) as context:
            taxon_graph.find_taxIDs(9)
        self.assertTrue("The parent_child_graph is empty. Check if file is correct and if first function "
                        "'create_graph' is called." in str(context.exception))

    def test_find_renamed_taxID(self):
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        self.assertEqual({2, 5, 6, 14, 15, 16, 37, 38}, taxon_graph.find_taxIDs(100))

    def test_find_non_existent_taxID(self):
        final_taxIDs = {6, 14, 15, 16}
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        final_taxIDs.update(taxon_graph.find_taxIDs(1000))
        self.assertEqual({6, 14, 15, 16}, final_taxIDs)

        with self.assertLogs('TaxonGraph', level='INFO') as log:
            taxon_graph.find_taxIDs(1000)
            self.assertEqual([
                'ERROR:TaxonGraph:User given taxon ID 1000 does not exist and is excluded from further analysis.'],
                log.output)

    def test_find_level_up(self):
        # taxID 38 from level subspecies to species to family without genus
        taxIDs = {2, 22, 18, 38}
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        taxIDs = {taxon_graph.find_level_up(taxID, 'genus') for taxID in taxIDs}
        self.assertEqual(taxIDs, {2, 3, 14})

    def test_child_until_species(self):
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        taxIDs = taxon_graph.find_taxIDs(3, True)
        self.assertEqual(taxIDs, {3,7,8,9,17,18,19,20,21})

    def test_find_next_common_ancestor(self):
        taxon_graph = TaxonGraph()
        taxon_graph.create_graph('./testdata/taxdumptest.tar.gz')
        taxonID = taxon_graph.find_next_common_ancestor([34, 30])
        self.assertEqual(taxonID, 4)


if __name__ == '__main__':
    unittest.main()
