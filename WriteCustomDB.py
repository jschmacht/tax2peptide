import gzip
import re
import logging
from pathlib import Path
import multiprocessing as mp
import shutil
import os
from tqdm import tqdm

logger = logging.getLogger(__name__)


class WriteCustomDB:
    def __init__(self, path_to_file, path_custom_db, taxonIDs=None):
        """
        :param path_to_file: path to the peptide database
        :param path_custom_db: path to the new generated taxon specific peptide database in fasta format )
        """
        self.path_to_file = Path(path_to_file)
        self.path_custom_db = path_custom_db
        self.taxonIDs = taxonIDs
        self.with_taxon_ID = False

    # open database file
    def open_database(self, type):
        """
        :param type: opening type: 'w', 'wb', 'r', 'rb'...
        : return opened file
        """
        try:
            f = open(str(self.path_to_file), type)
        except FileNotFoundError:
            logger.exception("Path to database does not exist.", exc_info=True)
            exit(1)
        return f

    # divide db-file in blocks for count of line numbers, if multiprocessing not possible
    def blocks(self, file, size=65536):
        while True:
            b = file.read(size)
            if not b: break
            yield b

    # compare all accession IDs in header line with set of taxon ID matching accession IDs
    def lookup_method_ncbi(self, line):
        """
        :param line: tuple of numbers, chunk[0]=entry point, chunk[1]= size to be readed
        """
        if '\x01' in line:
            # multiple_accessionIDs = {re.match("(.*?) ", hdr).group(1) for hdr in line[1:].split('\x01')}
            multiple_accessionIDs = {hdr.split(' ')[0] for hdr in line[1:].split('\x01')}
            for accessionID in multiple_accessionIDs:
                if accessionID in accessionIDs:
                    return True
            return False
        else:
            accessionID = (line[1:].split(' ')[0])
            if accessionID in accessionIDs:
                return True
            else:
                return False

    # ncbi peptide database: reads ncbi-db and write matching fasta entrys to taxon specified database (path_custom_db)
    def read_ncbi_linear(self, database):
        """
        :param database: opened peptide_database(nr.gz)
        """
        id_in = False
        try:
            with open(str(self.path_custom_db), 'w') as out:
                for line in tqdm(iter(database.readline, ''), total=self.linenumber):
                    line = line.rstrip()
                    if not line.startswith('>'):
                        if id_in:
                            out.write(line + '\n')
                    else:
                        if self.lookup_method_ncbi(line):
                            id_in = True
                            out.write(line + '\n')
                        else:
                            id_in = False
        except UnicodeDecodeError:
            logger.exception(
                "Failed opening new generated path '%s' to user specified taxon database." % self.path_custom_db)
    
    # ncbi peptide database: reads ncbi-db and write matching fasta entrys to taxon specified database (path_custom_db), use multiprocessing            
    def read_ncbi_mp(self, chunk):
        """
        :param chunk: tuple of numbers, chunk[0]=entry point in database file, chunk[1]= size to be readed
        """
        f = self.open_database("rb")
        f.seek(chunk[0])
        id_in = False
        fastas = ''
        for line in f.read(chunk[1]).decode("utf-8").splitlines():
            line = line.rstrip()
            if not line.startswith('>'):
                if id_in:
                    fastas += (line + '\n')
            else:
                if self.lookup_method_ncbi(line):
                    id_in = True
                    fastas += (line + '\n')
                else:
                    id_in = False
        return fastas

    # uniprot/swiss/trembl peptide database: reads ndb and write matching fasta entrys to taxon specified database (path_custom_db), use multiprocessing
    def read_uniprot_mp(self, chunk):
        """
        :param chunk: tuple of numbers, chunk[0]=entry point in database file, chunk[1]= size to be readed
        """
        f = self.open_database("rb")
        f.seek(chunk[0])
        re_expr = re.compile(r' *OX=(\d+)')
        id_in = False
        fastas = ''
        for line in f.read(chunk[1]).decode("utf-8").splitlines():
            line = line.rstrip()
            if not line.startswith('>'):
                if id_in:
                    fastas += (line + '\n')
            else:
                if int(re.search(re_expr, line).groups(1)[0]) in self.taxonIDs:
                    id_in = True
                    fastas += (line + '\n')
                else:
                    id_in = False
        f.close()
        return fastas

    # chunks have to end exact at end of a fasta entry
    def read_in_chunks(self, size=1024 * 1024):
        f = self.open_database("r")
        while True:
            start = f.tell()
            f.seek(size + start)
            s = f.readline()
            end = (f.tell() - start)
            while s:
                s = f.readline()
                if s.startswith('>') or not s:
                    break
                end = (f.tell() - start)
            f.seek(start + end)
            yield start, end
            if not s:
                f.close()
                break

    # divide database file in chunks for parallel processing
    def divide_into_chunks(self):
        chunks = []
        logger.debug('Start dividing database file in chunks.')
        for chunk in self.read_in_chunks():
            chunks.append(chunk)
        logger.debug('Database divided into %d chunks.' % len(chunks))
        return chunks
      
    # read database and store position information in list (value) to taxon_ID
    # if database is without taxon IDs, before matching scientifc name to ID       
    def read_database(self, with_taxon_ID, gzipped=False, threads=None, accessions=None):
        """   
         :param with_taxon_ID: uniprot or ncbi (with or without taxID in header) database
         :param gzipped: True if database file is gzipped
         :param threads: number of threads for multiprocessing
         :param accessions: set of accessionIDs, which match searched taxon IDs (for ncbi db)
         """

        # if database is .gz format: first unzipping:
        if gzipped:
            logger.info('Temporarily unzipping %s.' % self.path_to_file)
            path_to_gz_file = self.path_to_file
            self.path_to_file = self.path_to_file.parents[0] / (str(self.path_to_file.stem) + '.tmp')
            with gzip.open(str(path_to_gz_file), 'rb') as f_in:
                with open(str(self.path_to_file), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.debug('Unzipping done.')

        if not threads:
            threads = mp.cpu_count()

        # uniprot db
        if with_taxon_ID:
            chunks = self.divide_into_chunks()
            pool = mp.Pool(threads)
            try:
                with open(str(self.path_custom_db), 'w') as out:
                    # here multiprocessing starts
                    i, j = 0, 0
                    ten = int(len(chunks) / 10)
                    for fasta in pool.imap_unordered(self.read_uniprot_mp, chunks):
                        if i == ten:
                            j += 1
                            print('%d0%% readed.' % j)
                            i = 0
                        out.write(fasta)
                        i += 1
                pool.close()
                pool.join()
            except FileNotFoundError:
                logger.exception("Path to database (%s) does not exist." %str(self.path_custom_db), exc_info=True)
                exit(1)
        # ncbi db        
        else:
            global accessionIDs
            accessionIDs = accessions
            # in windows no forking possible and shared accessionIDs are to big to copy to different processes, in win only linear possible
            if not os.name == 'posix':
                try:
                    logger.debug('Count lines...')
                    with open(str(self.path_to_file), "rt") as database:
                        self.linenumber = sum(bl.count("\n") for bl in self.blocks(database))
                    logger.info('Reading/Writing...')
                    with open(str(self.path_to_file), "rt") as database:
                            self.read_ncbi_linear(database)
                except FileNotFoundError:
                    logger.exception("Path to database (%s) does not exist." %str(self.path_custom_db), exc_info=True)
                    exit(1)
          # ncbi, linux, multiprocessing
            else:
                chunks = self.divide_into_chunks()
                pool = mp.Pool(threads)
                try: 
                    with open(str(self.path_custom_db), 'w') as out:
                        # here multiprocessing starts
                        i, j = 0, 0
                        ten = int(len(chunks) / 10)
                        for fasta in pool.imap_unordered(self.read_ncbi_mp, chunks):
                            if i == ten:
                                j += 1
                                print('%d0%% readed.' % j)
                                i = 0
                            out.write(fasta)
                            i += 1
                    pool.close()
                    pool.join()
                except FileNotFoundError:
                    logger.exception("Path to database (%s) does not exist." % str(self.path_custom_db), exc_info=True)
                    exit(1)
        
        # remove tmp dezipped file
        if gzipped:
            self.path_to_file.unlink()


