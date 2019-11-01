import tarfile
import logging.config
from collections import defaultdict

logger = logging.getLogger(__name__)

class TaxonGraph():

    def __init__(self):
        self.parent_child_graph = defaultdict(list)
        self.child_rank_dict = {}
        self.oldtax_newtax_dict = {}
        self.taxon_name_dict = {}
        self.scientific_alternative_name_dict = {}
        self.child_parent_dict = {}
        self.synonym_taxon_dict = {}
        self.order = {'no rank': 0, 'subspecies': 1, 'varietas': 1, 'forma': 1, 'species': 2, 'species subgroup': 3,
                 'species group': 4, 'series': 5, 'subsection': 6, 'section': 7, 'subgenus': 8, 'genus': 9,
                 'subtribe': 10, 'tribe': 11, 'subfamily': 12, 'family': 13, 'superfamily': 14, 'parvorder': 15,
                 'infraorder': 16, 'suborder': 17, 'order': 18, 'superorder': 19, 'subcohort': 20, 'cohort': 21,
                 'infraclass': 22, 'subclass': 23, 'class': 24, 'superclass': 25, 'subphylum': 26, 'phylum': 27,
                 'superphylum': 28, 'subkingdom': 29, 'kingdom': 30, 'superkingdom': 31}

    # create taxon graph from nodes.dmp (in taxdump.tar.gz archive), a dictionary with taxon-ID-species name 
    # and a dictionary with changed/merged taxon IDs
    def create_graph(self, path_to_taxdump):
        """
        :param path_to_taxdump: path to location of the taxdump.tar.gz file (contains information of ancestry tree, source:NCBI)
          """
        # all lines split at | take 1. as parent, 0. as child and 3. as rank
        try:
            tar = tarfile.open(str(path_to_taxdump), "r:gz")
            nodes_handle = tar.extractfile("nodes.dmp")
            names_handle = tar.extractfile("names.dmp")
            merged_handle = tar.extractfile("merged.dmp")

            for line in nodes_handle:
                fields = [item.strip('\t') for item in line.decode(tarfile.ENCODING).split('|')]
                self.child_rank_dict[int(fields[0])] = fields[2]
                self.child_parent_dict[int(fields[0])] = int(fields[1])
            for line in merged_handle:
                fields = [item.strip('\t') for item in line.decode(tarfile.ENCODING).split('|')]
                self.oldtax_newtax_dict[int(fields[0])] = int(fields[1])
            for line in names_handle:
                fields = [item.strip('\t') for item in line.decode(tarfile.ENCODING).split('|')]
                if fields[3] == 'scientific name':
                    self.taxon_name_dict[int(fields[0])] = fields[1]
                else:
                    self.synonym_taxon_dict[fields[1]] = int(fields[0])

            for key, value in self.child_parent_dict.items():
                self.parent_child_graph[value].append(key)

            # delete circles in graph for top node
            for child, parent in self.child_parent_dict.items():
                if child == parent:
                    self.parent_child_graph[parent].remove(parent)
            # add leaves
            leaves = set(list(self.child_parent_dict.keys()) + list(self.child_parent_dict.values())) - set(self.parent_child_graph.keys())
            self.parent_child_graph.update(dict.fromkeys(leaves, []))

        except tarfile.ReadError as e:
            logger.exception("Can not read taxdump.tar.gz File. Without taxon tree, program quits.", exc_info=True)
            exit(1)

    # find for taxID all child tax IDs and returns them
    def find_taxIDs (self, taxID, childs_until_species=False):
        """
         :param taxID: tax ID for searching all children
         :type taxID:  positive int
         :return: all child taxon IDs of one taxon ID
         """
        if self.parent_child_graph == {}:
            raise Exception("The parent_child_graph is empty. Check if file is correct and if first function "
                            "\'create_graph\' is called.")
        species_taxIDs = set()
        # check if taxon IDs exist and have valid values
        if taxID not in self.child_rank_dict.keys():
            if self.oldtax_newtax_dict.get(taxID) is None:
                logger.error(
                    'User given taxon ID %d does not exist and is excluded from further analysis.' % taxID)
                return species_taxIDs
            taxID = self.oldtax_newtax_dict.get(taxID)

        if not childs_until_species:
            temp_children = {taxID}
            while temp_children:
                taxID = temp_children.pop()
                species_taxIDs.add(taxID)
                temp_children = temp_children.union(set(self.parent_child_graph.get(taxID)))
        else:
            # if rank of taxon ID is above species level, return same taxID, else go up until level 'species'
            taxID = self.find_level_up(taxID, 'species')
            temp_children = {taxID}
            while temp_children:
                taxID=temp_children.pop()
                if self.child_rank_dict[taxID] == 'species':
                    species_taxIDs.add(taxID)
                else:
                    species_taxIDs.add(taxID)
                    temp_children = temp_children.union(set(self.parent_child_graph.get(taxID)))

        return species_taxIDs

    # find for every taxID the given level up, if level not exist (this is regular) take the next smaller one
    # specified in order. If level of given taxID is higher as searched level, give taxID back unchanged
    def find_level_up(self, taxID, level):
        """
            :param taxID: tax ID for searching all children , positive int
            :param level: user given level for taxon search: string, one out of the strings in order dict
            :return taxon ID of specified level
                 """
        #check if level is really up, else take original taxID
        taxIDnew = taxID
        while self.child_rank_dict[taxIDnew] == 'no rank':
            taxIDnew = self.child_parent_dict.get(taxIDnew)
        # start taxID level higher, returned unchanged taxID
        if self.order[self.child_rank_dict[taxIDnew]] > self.order[level]:
            return taxID
        # start taxID level equals
        elif self.order[self.child_rank_dict[taxIDnew]] == self.order[level]:
            return taxIDnew
        # start taxID level is lower as wished
        else:
            last_specified_level_taxID = None
            while self.order[self.child_rank_dict[taxIDnew]] < self.order[level]:
                if self.child_rank_dict[taxIDnew] != 'no rank':
                    last_specified_level_taxID = (self.child_rank_dict[taxIDnew], taxIDnew)
                taxIDnew = self.child_parent_dict.get(taxIDnew)
            if self.order[self.child_rank_dict[taxIDnew]] == self.order[level]:
                return taxIDnew

            # sometimes not all pylogenetic ranks are filled (101570(=species) -> 196894(=no rank) -> 10699(=family))
            # rank genus and subfamily not present, take next highest with unique classification
            if self.order[self.child_rank_dict[taxIDnew]] > self.order[level]:
                if last_specified_level_taxID:
                    logger.warning("No ancestor with level %s of taxID %d exists. TaxID %d of level %s is "
                                   "returned." % (level, taxID, last_specified_level_taxID[1],
                                                  self.child_rank_dict[last_specified_level_taxID[1]]))
                    return last_specified_level_taxID[1]
                else:
                    logger.warning("Taxon ID %d has higher level than specified. Unchanged taxon ID %d of level %s is "
                                   "returned." % (taxID, taxID, self.child_rank_dict[taxID]))
                    return taxID

    # finds the next common ancestror of a set of taxon IDs
    def find_next_common_ancestor(self, taxonIDs):
        """
        :param taxonIDs: set or list of taxon IDs
        :return common ancestor of all taxon IDs
                 """
        # 1. find rank
        ranks=[]
        for taxon in taxonIDs:
            taxIDnew = taxon
            while self.child_rank_dict[taxIDnew] == 'no rank':
                taxIDnew = self.child_parent_dict.get(taxIDnew)
            ranks.append(self.order[self.child_rank_dict[taxIDnew]])
        index = ranks.index(max(ranks))
        highest_taxon = taxonIDs[index]
        del ranks[index]
        del taxonIDs[index]

        #list all parents:
        ancestor_list = [highest_taxon]
        while highest_taxon != 1:
            highest_taxon = self.child_parent_dict[highest_taxon]
            ancestor_list.append(highest_taxon)

        while taxonIDs:
            taxon = taxonIDs.pop()
            while taxon not in ancestor_list:
                taxon = self.child_parent_dict[taxon]
            # delete all list elem before found common ancestor
            index = ancestor_list.index(taxon)
            ancestor_list=ancestor_list[index:]

        return taxon