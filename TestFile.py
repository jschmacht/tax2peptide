import gzip
import logging.config
import re

logger = logging.getLogger(__name__)

class TestFile:

    # tests of taxon ID is in first header first line of database file and if file is gzipped
    # test taxon ID with reg exp ' OX='
    @staticmethod
    def test_gzipped(path_to_file):

        try:
            with open(path_to_file, 'rb') as f:
                b = f.read(2)
            if b != b'\037\213':
                gzipped = False
            else:
                gzipped = True
        except FileNotFoundError:
            logger.exception("Path to database does not exist.", exc_info=True)
            exit(1)
        try:
            with gzip.open(path_to_file, "rt") as db:
                db.readline()
        except OSError:
            gzipped = False
        return gzipped

    @staticmethod
    def test_uniprot(path_to_file):
        if TestFile.test_gzipped(path_to_file):
            try:
                with gzip.open(path_to_file, "rt") as db:
                    if re.search(' OX=', db.readline()):
                        with_taxon_ID = True
                    else:
                        with_taxon_ID = False
                return with_taxon_ID
            except FileNotFoundError:
                logger.exception("Path to database does not exist.", exc_info=True)
                exit(1)
        else:
            try:
                with open(str(path_to_file), "rt") as db:
                    if re.search(' OX=', db.readline()):
                        with_taxon_ID = True
                    else:
                        with_taxon_ID = False
                return with_taxon_ID
            except FileNotFoundError:
                logger.exception("Path to database does not exist.", exc_info=True)
                exit(1)
