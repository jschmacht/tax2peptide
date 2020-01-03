import logging.config
import gzip
import multiprocessing as mp
import shutil
from pathlib import Path
from TestFile import TestFile

logger = logging.getLogger(__name__)

class Accession:

    def __init__(self, taxIDs):
        """
            :param taxIDs: set of all taxIDs to search with
        """
        self.taxIDs = taxIDs
        self.accessionIDs = set()
        self.prot_path = None

    # return generators, chunks end at end of complete lines
    def read_in_chunks(self, file, size=1024 * 1024):
        """
            :param file: path to database
            :param size: size of chunks
            :return: generators of entry point and read length
        """
        f = open(str(file), 'r')
        f.readline()
        while True:
            start = f.tell()
            f.seek(size + start)
            s = f.readline()
            end = (f.tell() - start)
            f.seek(start + end)
            yield start, end
            if not s:
                f.close()
                break

    # worker function of multiprocessing part, read part of accession file and insert matching accession IDs in part_list
    def read_chunks_into_part_list(self, chunk):
        """
                :param chunk: tuple of numbers, chunk[0]=entry point, chunk[1]= size to be readed
                :return: list with accession IDs
                  """
        try:
            f = open(str(self.prot_path), "rb")
            f.seek(chunk[0])
            part_list = []
            for line in f.read(chunk[1]).decode("utf-8").splitlines():
                fields = [item.strip('\t') for item in line.split('\t')]
                if int(fields[2]) in self.taxIDs:
                    part_list.append(fields[1])
            f.close()
            return part_list
        except FileNotFoundError:
            logger.exception("Path to database does not exist.", exc_info=True)
            exit(1)

    # read accessions, if taxon ID matches to taxon set, accession added to self.accessionIDs set.
    def read_accessions(self, protaccession_path, pdbaccession_path, threads=None):
        """
        :param protaccession_path: path to location of the prot.accession2taxid.gz file
        (contains information protein accession numbers and their taxon IDs, source:NCBI)
        :param pdbaccession_path: path to location of the pdb.accession2taxid.gz file
        (contains information pdb accession numbers and their taxon IDs, source:NCBI)
        :param threads: number of threads to use for multiprocessing
          """
        # if prot.accession2taxid.gz gzipped -> first unzip
        gzipped = TestFile.test_gzipped(protaccession_path)
        if gzipped:
            logger.info('Temporarily unzipping prot2accession.gz file.')
            self.prot_path = Path(protaccession_path).parents[0] / (str(Path(protaccession_path).stem) + '.tmp')
            with gzip.open(str(protaccession_path), 'rb') as f_in:
                with open(str(self.prot_path), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info('Unzipping done.')
        else:
            self.prot_path = protaccession_path

        if not threads:
            threads = mp.cpu_count()
        pool = mp.Pool(threads)
        # chunks = list of all chunks list tuples: (entry point to db, number of chars to be read)
        chunks = []
        for chunk in self.read_in_chunks(self.prot_path):
            chunks.append(chunk)

        logger.info('Start reading accession2prot database file with %s threads.' % str(threads))
        # here multiprocessing starts, results saved in self.accessionIDs
        accession_list = []
        i, j = 0, 0
        ten = int(len(chunks)/10)
        for part_list in pool.imap(self.read_chunks_into_part_list, chunks):
            if i == ten:
                j += 1
                print('%d0%% readed.' % j)
                i = -1
            accession_list.append(part_list)
            i += 1
        pool.close()
        for l in accession_list:
            self.accessionIDs.update(set(l))

        # delete tmp unzipped prot.accession2taxid file
        if gzipped:
            self.prot_path.unlink()
        # reading and procession small pdbaccession file, no multiprocession
        try:
            with gzip.open(str(pdbaccession_path), "rt") as accession_handle:
                accession_handle.readline()
                for line in accession_handle:
                    fields = [item.strip('\t') for item in line.split('\t')]
                    if int(fields[2]) in self.taxIDs:
                        self.accessionIDs.add(fields[1])
            logger.debug('Number of matching accession IDs: %d' % len(self.accessionIDs))
        except FileNotFoundError:
            logger.exception("Path to accession file %s does not exist." % pdbaccession_path, exc_info=True)
            exit(1)
        except OSError:
            logger.exception("File %s has wrong format. Need gzipped file. " % protaccession_path, exc_info=True)
            exit(1)
