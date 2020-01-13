"""
#
#-----------------------------
#  Robert Koch Institut
#         2019
#-----------------------------
# Author:
#
# Juliane Schmachtenberg
#
#------------------------------
"""

import re
import argparse
import urllib.request
from datetime import date
import logging
from Download import Download
from TaxonGraph import TaxonGraph
from DatabaseCleaner import DatabaseCleaner
from Accession import Accession
from TestFile import TestFile
from WriteCustomDB import WriteCustomDB
from Output import Output
import shutil
import pickle
from pathlib import Path
from version import __version__


# Initializing logger, log file saved in direction output_path = same as database
def initialize_logger(output_path, verbose=None):
    logpath = output_path.parents[0]
    logfile_name = str(output_path.stem + '.log')
    try:
        Path(logpath / logfile_name).unlink()
    except FileNotFoundError:
        pass
    if verbose:
        logging.basicConfig(format='%(asctime)s | %(name)s | %(levelname)s: %(message)s',
                        level=logging.DEBUG,
                        handlers=[logging.FileHandler("{0}/{1}".format(logpath, logfile_name)), logging.StreamHandler()])
    else:
        h1 = logging.FileHandler("{0}/{1}".format(logpath, logfile_name))
        h1.setLevel(logging.DEBUG)
        h2 = logging.StreamHandler()
        h2.setLevel(logging.INFO)
        logging.basicConfig(format='%(asctime)s | %(name)s | %(levelname)s: %(message)s',
                            level=logging.DEBUG,
                            handlers=[h1, h2])

    return logging.getLogger(__name__)

# reads metadata to download files from uniprot, like md5-hash and db-version
def read_uniprot_metadata(url, database, logger):
    attempts = 0
    while attempts < 5:
        try:
            with urllib.request.urlopen(url) as f:
                metavalues_db = f.read().decode('ascii')
            hit_version = re.search("<version>(.*)</version>", metavalues_db)
            hit_hash = re.search('<file name="' + database + '">(.*\n){3} *<hash type="md5">(.*)</hash>',
                                 metavalues_db)
            version = hit_version.group(1) if hit_version else None
            md5 = hit_hash.group(2) if hit_hash else None
            return version, md5
        except urllib.error.URLError:
            attempts += 1
            if attempts == 5:
                logger.exception(
                    'Download of the uniprot database version and md5 hash for validation of the file %s failed. \n'
                    'Program continues without validation.' % database, exc_info=True)
                return None, None

# reads md5 hash of ncbi non redundant database nr.gz
def read_ncbi_hash(url, logger):
    attempts = 0
    while attempts < 5:
        try:
            with urllib.request.urlopen(url) as f:
                return f.read().split()[0].decode('ascii')
        except urllib.error.URLError:
            attempts += 1
            if attempts == 5:
                logger.exception(
                    'Download md5 hash from %s for validation of the ncbi file failed.Program continues '
                    'without validation.' % url, exc_info=True)
                return None

# check if commandline input is positive integer
def positive_integer(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value." % value)
    return ivalue


def main():
    parser = argparse.ArgumentParser(description='Find all protein database entrys of specified taxon IDs and their descendants.' \
          ' One taxID or a taxID input file must be provided. Peptide-Databases from NCBI or Uniprot can be used. User defined databases,' \
          ' if header contain taxon IDs (e.g. OX=1111) or ncbi/uniprot accession IDs.')
    parser.add_argument('-i', '--input', dest='input', default=None, help='TaxID input file: tabular file containing a column of NCBI'
                                                                          ' taxon IDs. Columns tab separated.')
    parser.add_argument('-c', '--column', dest='column', type=positive_integer, default=0, help='The column (zero-based) in the tabular '
                                                                                  'file that contains Taxon IDs. Default = 0.')
    parser.add_argument('-t', '--taxon', dest='taxon', type=positive_integer, nargs='+', action='append',
                        help='NCBI taxon ID/s for database extraction. Multiple taxonIDs seperated by space.')
    parser.add_argument('-d', '--database', dest='database', choices=['ncbi', 'uniprot', 'swissprot', 'trembl'],
                        default='uniprot', help='Database choice for analysis or for download. Choices: ncbi, uniprot, tremble, swissprot. '
                        'No download, if databases with original name are stored in same folder as option --path ')
    parser.add_argument('-p', '--path', dest='path', default=None, help='Path to folder with all needed '
                            'databases: taxdump.tar.gz (for all databases), prot.accession2taxid or prot.accession2taxid.gz and '
                            'pdb.accession2taxid.gz (for ncbi databases). Optional: peptide_database named: nr/nr.gz, '
                            'uniprot_trembl.fasta/uniprot_trembl.fasta.gz or uniprot_sprot.fasta/uniprot_sprot.fasta.gz'
                            ' or uniprot.fasta./uniprot.fasta.gz')
    parser.add_argument('-o', '--out', dest='out', default=None,
                        help="File name and direction of the result taxon specified peptide database. "
                             "Default = /taxon_specified_db_DATE/taxon_database.fasta")
    parser.add_argument('-n', '--dbname', dest='dbname', default=None,
                        help="Database name and direction. If database is in other folder than --path or name deviates from standard names.")
    parser.add_argument('-l', '--level', dest='level',  choices=['species', 'section', 'genus', 'tribe', 'subfamily', 'family', 'superfamily',
                                                                 'order', 'superorder', 'class', 'phylum', 'kingdom', 'superkingdom'], default=None,
                        help='Hierarchy level up in anchestral tree. Choices: species, section, genus, tribe, '
                             'subfamily, family, superfamily, order, superorder, class, phylum, kingdom, superkingdom')
    parser.add_argument('-z', '--no_descendants', dest='no_descendants', action='store_true', default=False,
                        help='Select peptide database only by given taxon IDs, descendant taxons are excluded.')
    parser.add_argument('-s', '--species', dest='species', action='store_true', default=False,
                        help='Select peptide database only until taxonomic level "species", descendents from species are excluded.')
    parser.add_argument('-r', '--non_redundant', dest='non_redundant', action='store_true', default=False,
                        help='Makes the final database non redundant in regard to sequences, headers are concatenated.')
    parser.add_argument('-u', '--threads', dest='threads', type=positive_integer, action="store",
                        help='Number of threads for using multiprocessing. Default = number of cores.')
    parser.add_argument('-x', '--reduce_header', dest='reduce_header', action='store_true', default=False,
                        help='Reduce the long headers of NCBI entries to accession IDs. Use only for NCBI databases.')
    parser.add_argument('--version', action='version', version=('version ' + __version__))
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                        help='Verbose shows details about program progress and more information.')

    options = parser.parse_args()
    # url adresses for download:
    url_protaccession2taxID = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz'
    url_protaccession2taxID_md5 = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5'
    url_pdbaccession2taxID = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz'
    url_pdbaccession2taxID_md5 = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz.md5'
    url_taxdump = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    url_taxdump_md5 = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5'
    url_database_ncbi = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz'
    url_database_md5_ncbi = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5'
    url_database_swissprot = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
    url_database_trembl = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz'
    url_uniprot_metadata = 'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/RELEASE.metalink'
    db_dict_name = {'ncbi': url_database_ncbi.split('/')[-1], 'uniprot': 'uniprot.fasta.gz',
               'swissprot': url_database_swissprot.split('/')[-1], 'trembl': url_database_trembl.split('/')[-1]}

    # if not option out, a new folder with name taxon_database and date for result database and log file is created
    if options.out:
        output_path = Path.cwd() / options.out
    else:
        output_path = Output.createDir(Path.cwd())

    logger = initialize_logger(output_path, options.verbose)

    # if path to database given, check uniprot or ncbi format
    if options.dbname:
        if not TestFile.test_uniprot(options.dbname):
            options.database = 'ncbi'

    for arg, value in sorted(vars(options).items()):
        logger.debug("Argument %s: %r", arg, value)
    logger.info("Result database and log file are saved in direction %s" % output_path)

    # try open config file and read path to database folder, if no path option is entered
    config = False
    if not options.path:
        try:
            path_to_main =Path(__file__, '..').resolve()
            with open(str(path_to_main) + "/tax2peptide.config", 'r') as config:
                database_folder = Path(config.readline().strip())
                config = True
            path_to_db = database_folder/db_dict_name[options.database]
            taxdump_b, prot_gz_b, prot_b, pdb_b, db_gz_b, db_b = Output.check_files_exist(
                [database_folder / url_taxdump.split('/')[-1], database_folder / url_protaccession2taxID.split('/')[-1],
                  database_folder/'prot.accession2taxid', database_folder / url_pdbaccession2taxID.split('/')[-1],
                 path_to_db, path_to_db.parents[0]/path_to_db.stem])
            if db_b: path_to_db = path_to_db.parents[0] / path_to_db.stem
            if not taxdump_b:
                logger.warning("File taxdump.tar.gz does not exist does not exist under the path %s"
                               " and will be downloaded." % str(database_folder))
            if not pdb_b and options.database == 'ncbi':
                logger.warning("File pdb.accession2taxid.gz does not exist does not exist under the path %s"
                               " and will be downloaded." % str(database_folder))
            if not prot_gz_b and not prot_b and options.database == 'ncbi':
                logger.warning("File prot.accession2taxid.gz does not exist does not exist under the path %s"
                               " and will be downloaded." % str(database_folder))
            if options.dbname is None and not db_b and not db_gz_b:
                logger.warning("Database file %s does not exist does not exist under the path %s"
                               " and will be downloaded." % (db_dict_name[options.database], str(database_folder)))
        except FileNotFoundError:
            pass
    # if options.path specified: folder to all databases (without peptide db if options.dbname)
    # check if all needed files in database folder: bool values _b: True = file exists and not downloaded again
    # if not exist , create folder with user defined name in option --path
    else:
        database_folder = Path.cwd() / options.path
        if database_folder.exists():
            path_to_db = database_folder / db_dict_name[options.database]
            taxdump_b, prot_gz_b, prot_b, pdb_b, db_gz_b, db_b = Output.check_files_exist([database_folder/url_taxdump.split('/')[-1],
                database_folder/url_protaccession2taxID.split('/')[-1], database_folder/'prot.accession2taxid',
                database_folder/url_pdbaccession2taxID.split('/')[-1], path_to_db, path_to_db.parents[0]/path_to_db.stem])
            if db_b:
                path_to_db = path_to_db.parents[0]/path_to_db.stem
            if not taxdump_b:
                logger.warning("File taxdump.tar.gz does not exist does not exist under the path %s"
                               " and will be downloaded." % str(database_folder))
            if not pdb_b:
                logger.warning("File pdb.accession2taxid.gz does not exist does not exist under the path %s"
                               " and will be downloaded." % str(database_folder))
            if not prot_gz_b and not prot_b:
                logger.warning("File prot.accession2taxid.gz does not exist does not exist under the path %s"
                               " and will be downloaded." % str(database_folder))
            if options.dbname is None and not db_b and not db_gz_b:
                logger.warning("Database file %s does not exist does not exist under the path %s"
                               " and will be downloaded." % (db_dict_name[options.database], str(database_folder)))
        else:
            try:
                database_folder.mkdir()
                logger.exception("Folder of option --path (%s) does not exist and is now new generated. All needed files "
                                 "are downloaded and stored in given direction. " % database_folder, exc_info=True)
                prot_gz_b = prot_b = pdb_b = taxdump_b = db_gz_b = db_b = False
            except OSError:
                logger.exception("Folder of option --path (%s) does not exist and con not be created." % database_folder,
                                 exc_info=True)
                exit(1)

    # if options.path not given, new folder for all needed databases, check if new generated path exists, check content folder
    if options.path is None and not config:
        database_folder = Path.cwd() / ('databases_' + str(date.today()))
        try:
            database_folder.mkdir()
            logger.info("Downloaded databases are saved in direction %s" % database_folder)
            taxdump_b = prot_gz_b = prot_b = pdb_b = db_gz_b = db_b = False
        except FileExistsError:
            path_to_db = database_folder/db_dict_name[options.database]
            logger.debug("Database folder %s already exists. Checking for content." % database_folder)
            if options.database == 'ncbi':
                taxdump_b, prot_gz_b, prot_b, pdb_b, db_gz_b, db_b = Output.check_content([database_folder/url_taxdump.split('/')[-1],
                database_folder/url_protaccession2taxID.split('/')[-1], database_folder/'prot.accession2taxid',
                database_folder/url_pdbaccession2taxID.split('/')[-1], path_to_db, path_to_db.parents[0]/path_to_db.stem])
                if db_b: path_to_db = path_to_db.parents[0] / path_to_db.stem
            else:
                taxdump_b, db_gz_b, db_b = Output.check_content([database_folder / url_taxdump.split('/')[-1], path_to_db, path_to_db.parents[0]/path_to_db.stem])
                prot_gz_b = prot_b = pdb_b = False
                if db_b: path_to_db = path_to_db.parents[0] / path_to_db.stem
        except OSError:
            logger.exception("No permission to create new folders.", exc_info=True)
            exit(1)


    # given path to database checked, if not exists quit.
    if options.dbname:
        path_to_db = Path.cwd()/options.dbname
        db_b = Output.check_files_exist([path_to_db])[0]
        if not db_b:
            logger.error('Given database %s does not exist. Enter correct path under option --dbname. Program quits.' % path_to_db)
            exit(1)

    # download all databases same day:
    # download taxdump file
    if not taxdump_b:
        taxdump_md5 = read_ncbi_hash(url_taxdump_md5, logger)
        dwl_taxdb = Download(url_taxdump, database_folder/url_taxdump.split('/')[-1], taxdump_md5)
        dwl_taxdb.download()
        logger.debug('End download of taxdump.tar.gz')
    # download prot.accession2taxid.gz (only for ncbi) and check md5 hash
    if not prot_gz_b and not prot_b and options.database == 'ncbi':
        md5_hash = read_ncbi_hash(url_protaccession2taxID_md5, logger)
        dwl_protaccession = Download(url_protaccession2taxID, database_folder/url_protaccession2taxID.split('/')[-1],
                                     md5=md5_hash)
        dwl_protaccession.download()
        logger.debug('End download from %s to location %s.' %
                    (url_protaccession2taxID, str(database_folder/url_protaccession2taxID.split('/')[-1])))
    # download pdb.accession2taxid.gz (only for ncbi) and check md5 hash
    if not pdb_b and options.database == 'ncbi':
        md5_hash = read_ncbi_hash(url_pdbaccession2taxID_md5, logger)
        dwl_pdbaccession = Download(url_pdbaccession2taxID, database_folder / url_pdbaccession2taxID.split('/')[-1],
                                    md5=md5_hash)
        dwl_pdbaccession.download()
        logger.info('End download from %s to location %s.'
                    % (url_pdbaccession2taxID, str(database_folder / url_pdbaccession2taxID.split('/')[-1])))
    # download peptide database and check md5 hash
    if not db_b and not db_gz_b:
        if options.database == 'ncbi':
            database_version_ncbi = 'ncbi ' + str(date)
            md5_hash = read_ncbi_hash(url_database_md5_ncbi, logger)
            dwl_db = Download(url_database_ncbi, database_folder / db_dict_name['ncbi'], md5=md5_hash)
            dwl_db.download()
            logger.info("Databaseversion: %s" % database_version_ncbi)
            path_to_db = database_folder / db_dict_name['ncbi']
        else:
            if options.database == 'swissprot' or options.database == 'uniprot':
                database_version_swissprot, hash_swissprot = read_uniprot_metadata(url_uniprot_metadata, db_dict_name['swissprot'],logger)
                logger.info("Database version swissprot: %s " % database_version_swissprot)
                dwl_db_swiss = Download(url_database_swissprot, database_folder / db_dict_name['swissprot'], md5=hash_swissprot)
                dwl_db_swiss.download()
                path_to_db = database_folder / db_dict_name['swissprot']
            if options.database == 'trembl' or options.database == 'uniprot':
                database_version_trembl, hash_trembl = read_uniprot_metadata(url_uniprot_metadata, db_dict_name['trembl'],logger)
                logger.info("Databaseversion trembl: %s." % database_version_trembl)
                dwl_db_trembl = Download(url_database_trembl, database_folder / db_dict_name['trembl'], md5=hash_trembl)
                dwl_db_trembl.download()
                path_to_db = database_folder / db_dict_name['trembl']
            # concetenate  swissprot and trembl to uniprot file
            if options.database == 'uniprot':
                try:
                    logger.info("Concatenate swissprot and trembl to uniprot database with name uniprot.fasta")
                    with open(str(database_folder /db_dict_name['trembl']), 'ab') as trembl:
                        with open(str(database_folder /db_dict_name['swissprot']), 'rb') as swissprot:
                            shutil.copyfileobj(swissprot, trembl)
                    # rename trembl to uniprot:
                    Path(database_folder / db_dict_name['trembl']).rename(database_folder / db_dict_name['uniprot'])
                    logger.info("Uniprot database is now ready.")
                    path_to_db = database_folder / db_dict_name['uniprot']
                except FileNotFoundError:
                    logger.exception("Creation of uniprot database file out of trembl and swissprot file failed.",
                                     exc_info=True)
                    exit(1)

    # create config file
    try:
        path_to_main = Path(__file__, '..').resolve()
        with open(str(path_to_main / "tax2peptide.config"), 'w') as config:
            config.write(str(database_folder) + '\n')
    except OSError:
        logger.debug('Can not create config file')

    # Read taxIDs from option -t and option -i
    if options.taxon:
        taxIDs = set([taxID for taxonlist in options.taxon for taxID in taxonlist])
    else:
        taxIDs = set()
    if options.input:
        try:
            with open(options.input, 'r') as inputFile:
                for i, line in enumerate(inputFile):
                    fields = line.rstrip('\r\n').split('\t')
                    if len(fields) >= abs(options.column):
                        taxID = fields[options.column].strip()
                        if taxID.isdigit():
                            taxIDs.add(int(taxID))
                        else:
                            logger.error('Value %s in line %i of taxon input file is not a number. '
                                         'Right column number specified?' % (taxID, i))
                            continue
                    else:
                        logger.error('Column number is bigger as number of columns in taxon ID input file. '
                                     'Program continues without taxon IDs from input file.')
        except FileNotFoundError:
            logger.exception('Taxon ID input file does not exist under specified path.', exc_info=True)

    if not taxIDs:
        logger.error('No taxon ID given. Please check your input. Program quits. ')
        raise Exception('No taxon IDs.')
        exit(1)

    logger.debug('Given Tax-IDs: %s' % ' '.join(str(it) for it in taxIDs))

    # Try load pre-builded taxonomy graph or built taxonomy graph now
    if not (database_folder / 'taxon_graph').is_file():
        taxon_graph = TaxonGraph()
        logger.debug("Start building taxon graph.")
        taxon_graph.create_graph(database_folder/url_taxdump.split('/')[-1])
        logger.debug("Taxon graph successfully build.")
        # save TaxonGraph to harddrive:
        with open(str(database_folder / 'taxon_graph'), 'wb') as handle:
            pickle.dump(taxon_graph, handle, protocol=pickle.HIGHEST_PROTOCOL)
            logger.debug('Safe taxon graph to location: %s' % str(database_folder / 'taxon_graph'))
    # load Taxon Graph
    else:
        try:
            logger.debug('Load taxon graph.')
            with open(str(database_folder / 'taxon_graph'), 'rb') as handle:
                taxon_graph = pickle.load(handle)
        except UnicodeDecodeError or EOFError:
            logger.exception(
                "Failed opening path to taxon graph / taxon_graph is corrupted. Delete %s file."
                % str(database_folder / 'taxon_graph'))
            exit(1)

    # adjusts the hierarchy level, if level does not exist, take next smaller level
    if options.level:
        logger.debug("Start selection of next ancestor of level %s for all given taxIDs" % options.level)
        taxIDs = {taxon_graph.find_level_up(taxID, options.level) for taxID in taxIDs}
        logger.info("All taxon IDs are set up to level %s in anchestral tree. Taxon IDs of level %s: %s"
                    % (options.level, options.level, ' '.join(str(it) for it in taxIDs)))

    final_taxIDs = set()
    # find all descendants
    if not options.no_descendants:
        logger.debug("Start searching for all child taxon IDs.")
        for taxID in taxIDs:
            final_taxIDs.update(taxon_graph.find_taxIDs(taxID, options.species))
        logger.debug("End searching for all child taxon IDs.")
        logger.debug('Number of final taxon IDs: %s' % str(len(final_taxIDs)))
    else:
        final_taxIDs = taxIDs
        logger.debug('Number of taxon IDs for database search: %s' % str(len(final_taxIDs)))

    # generate accession_taxID dict for ncbi db search and write custom specified db to --out
    with_taxon_ID = TestFile.test_uniprot(path_to_db)
    if not with_taxon_ID:
        accession = Accession(final_taxIDs)
        logger.debug('Read accession files.')
        if prot_b:
            accession.read_accessions(database_folder/'prot.accession2taxid',
                                  database_folder/url_pdbaccession2taxID.split('/')[-1], options.threads)
        else:
            accession.read_accessions(database_folder/url_protaccession2taxID.split('/')[-1],
                                  database_folder/url_pdbaccession2taxID.split('/')[-1], options.threads)
        logger.debug('All accession IDs collected.')
        logger.info('Start writing taxon selected peptide database to %s.' % output_path)
        wc = WriteCustomDB(path_to_db, output_path)
        wc.read_database(False, gzipped=TestFile.test_gzipped(path_to_db), accessions=accession.accessionIDs,
                         threads=options.threads)
        logger.debug('End writing taxon selected peptide database.')
        # non redundant database

    # uniprot: write custom specified db to --out
    else:
        logger.info('Start writing taxon selected peptide database to %s.' % output_path)
        wc = WriteCustomDB(path_to_db, output_path, final_taxIDs)
        wc.read_database(True, threads=options.threads, gzipped=TestFile.test_gzipped(path_to_db))
        logger.debug('End writing taxon selected peptide database.')

    # non redundant database
    if options.non_redundant:
        DatabaseCleaner.non_redundant(output_path, with_taxon_ID)
        # remove redundant database:
        output_path.unlink()

    if options.reduce_header and not with_taxon_ID:
        # reduce headers of NCBI database
        DatabaseCleaner.reduce_header(output_path)
        output_path.unlink()


    logger.info('Program finished.')
    exit(0)


if __name__ == '__main__':
    main()
