# tax2peptide

tax2peptide creates based on the given taxon IDs and a reference database a taxon specific database in fasta format. This taxon specific database contains according to the selected options all fasta entries of:
* **the given taxon IDs and their descendant taxon IDs in the phylogenetic tree**
* **the given taxon IDs (option --no_descendants)** 
* **the given taxon IDs adapted to specified level up in the phylogenetic tree and their descendant taxon IDs (option --level)**
* **the given taxon IDs and their descendant taxon IDs in the phylogenetic tree until level species, taxonIDs with lower level are not included (option --species)**

Databases from which the matching entries are read out are the NCBI non-redundant peptide database, swissprot, uniprot or trembl database. Also user-defined databases are possible as long as the header of the database contains taxonIDs in form of "OX=NUMBER" or contains NCBI/uniprot accession numbers.
Using uncompressed databases would speed up the program considerably.

## Getting Started

Tax2Peptide is a python3 command line tool. It can be installed as pypi-package or as conda-package.

### Prerequisites

Python3

needed packages: tqdm, biopython, wget
```
pip install tqdm biopython wget
```

### Installing
```
pip install -i https://pypi.org/simple/tax2peptide-jschmacht 
```


## Deployment

Tax2Peptide is a command line tool and starts with: 
```
python3 path/to/tax2peptide -options
```

### Options:
-i 	--input   TaxID input file: tabular file containing a column of NCBI taxon IDs. Columns tab separated.

-c 	--column   The column (zero-based) in the tabular file that contains Taxon IDs. Default = 0.

-t 	--taxon   NCBI taxon ID/s for database extraction. Multiple taxonIDs seperated by space.

-d 	--database   Database choice for analysis or for download. Choices: ncbi, uniprot, tremble, swissprot.

-p 	--path   Path to folder with all needed databases: taxdump.tar.gz (for all databases), prot.accession2taxid or prot.accession2taxid.gz and pdb.accession2taxid.gz (for ncbi databases). Optional: peptide_database named: nr/nr.gz or uniprot_trembl.fasta/uniprot_trembl.fasta.gz or uniprot_sprot.fasta/uniprot_sprot.fasta.gz or uniprot.fasta./uniprot.fasta.gz

-o 	--out   File name and direction of the result taxon specified peptide database. Default = /taxon_specified_db_DATE/taxon_database.fasta

-n 	--dbname   Database name and direction. If database is in other folder than --path or have a different name.

-l 	--level   Hierarchy level up in anchestral tree. Choices: species, section, genus, tribe, subfamily, family, superfamily, order, superorder, class, phylum, kingdom, superkingdom

-r 	--non_redundant   Make the final database non redundant in regard to sequences, headers are concatenated.

-z 	--no_descendants   Select peptide database only by given taxon IDs, descendant taxons are excluded.

-s 	--species   Select peptide database only until taxonomic level "species", descendants from species are excluded.

-u 	--threads   Number of threads for using multiprocessing. Default = number of cores.

-x 	--reduce_header   Reduce the long headers of NCBI entries to accession IDs. Use only for NCBI databases.

### Dependencies:
Needed databases for NCBI reference database: 
* protaccession2tax.gz / protaccession2tax
* pdbaccession2tax.gz
* taxdump.tar.gz
* nr.gz / nr

Needed databases for uniprot/swissprot/trembl reference database: 
* taxdump.tar.gz
* uniprot.fasta.gz / uniprot.fasta / uniprot_sprot.fasta.gz / uniprot_sprot.fasta / uniprot_trembl.fasta.gz / uniprot_trembl.fasta

All database files should be downloaded the same day and stored in the same folder.

#### Databases
All should be downloaded at the same date as the peptide database to ensure successful accession matching.
The database can be downloaded manually or downloaded by Tax2Peptide with option --database {ncbi, uniprot, trembl, swissprot}, file integrity would be md5-checked.


| database name       | description                                                         |source | adress
| ------------------- |---------------------------------------------------------------------|-------| --------------------------------------------------------------------------------------------------------------------|
| NCBI                | non redundant peptide database                                      |NCBI   | ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz                                                                     |
| Swissprot           | curated peptide database                                            |Uniprot| ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz           |
| Trembl              | peptide database                                                    |Uniprot| ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz          |
| Uniprot             | concatenated swissprot and trembl database                          |Uniprot|                                                                                                                     |
| prot.accession2taxid|contain links between accession IDs and taxonomic lineage (taxon IDs)|NCBI   | ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz                                           |
| pdb.accession2taxid |contain links between accession IDs and taxonomic lineage (taxon IDs)|NCBI   | ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz                                            |
| taxdump             |tar-gz-compressed taxdump file containing information about the phylogenetic lineage and links between taxIDs and scientific names etc.|NCBI   |ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz            |

#### Usage with download: 
All needed databases will be downloaded to specified path (option --path). If not specified a folder with name databases_DATE will be used as default.

Examples of usage:
```
tax2peptide.py -d uniprot -i path/to/input/taxon_ID_file  -> new Folder databases_DATE with: taxdump.tar.gz, uniprot.fasta
```
```
tax2peptide.py -i path/to/input/taxon_ID_file  -> new Folder databases_DATE with: taxdump.tar.gz, uniprot.fasta
```
```
tax2peptide.py -d ncbi -p path/to/my_new_databases -i path/to/input/taxon_ID_file  -> new Folder/used Folder my_new_databases with: protaccession2tax.gz, pdbaccession2tax.gz, taxdump.tar.gz, nr.gz
```
#### Usage if all database files are already downloaded:
- option --path determines folder with all needed databases
(optional --dbname

- at least one taxon ID or taxon ID input file must be provided
Path will be checked for all required database files and missing databases will be downloaded
Examples of usage:
```
tax2peptide.py -p path/to/folder -n path/to/reference_database -t 11111 22222 -o path/my_taxon_specified_database.fasta
```
```
tax2peptide.py -p path/to/folder -n path/ to/ uniprot.fasta -t 11111 22222 -i path/to/input
```
```
tax2peptide.py -d ncbi -p path/to/folder -i path/to/input
```
```
tax2peptide.py -d uniprot -p path/to/folder -i path/to/input -o path/to/user_specified_db.fasta
```
If path is once determined, it must not be specified again, as long as the same folder shell be used.


## Authors

* **Juliane Schmachtenberg** 

[project_on_github](https://github.com/jschmacht/tax2peptide)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
