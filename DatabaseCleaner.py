import logging.config
from Bio import SeqIO
from pathlib import Path

logger = logging.getLogger(__name__)


class DatabaseCleaner:

    # modified https://biopython.org/wiki/Sequence_Cleaner
    @staticmethod
    def non_redundant(fasta_file, with_taxon_ID):
        """
        :param fasta_file: fasta file to make non redundant
        :param with_taxon_ID: True for uniprot, False for NCBI, different delimiter in header
        """
        logger.info('Start creating non redundant database.')
        sequences = {}
        nr_file = str(fasta_file.parents[0] / str(fasta_file.stem + '_nr.fasta'))
        # Using the Biopython fasta parse to read fasta input
        for seq_record in SeqIO.parse(str(fasta_file), "fasta"):
            sequence = str(seq_record.seq).upper()
            if sequence not in sequences:
                sequences[sequence] = seq_record.description
            # If sequence is already in the dict, concatenate headers to the other one that is already in the hash table
            else:
                if not with_taxon_ID:
                    sequences[sequence] += '\x01' + seq_record.description
                else:
                    sequences[sequence] += '|' + seq_record.description

        # Write non redundant sequences
        try:
            with open(str(nr_file), "w+") as output_file:
                # Just read the hash table and write on the file as a fasta format with line length 60
                for sequence in sequences:
                    # write header
                    output_file.write(">" + sequences[sequence] + "\n")
                    # write sequence
                    seq_parts = [sequence[i:i + 60] for i in range(0, len(sequence), 60)]
                    for seq in seq_parts:
                        output_file.write(seq + "\n")
        except OSError:
            logger.exception('Not able to write non redundant database.', exc_info=True)

        logger.info("Fasta database %s is now non redundant and saved to %s" % (fasta_file.name, nr_file))

    # Reduce headers to accessionIDs separated by \x01. file unzipped
    @staticmethod
    def reduce_header(path_to_db):
        """
                :param path_to_db: ncbi fasta file to reduce headers
        """
        try:
            with open(path_to_db, "rt") as database:
                rh_path = str(Path(path_to_db).parents[0] / str(Path(path_to_db).stem + '_rh.fasta'))
                with open(rh_path, 'w') as out:
                    logger.info('Start reducing NCBI headers.')
                    for line in database:
                        if line.startswith('>'):
                            if '\x01' in line:
                                out.write('>')
                                multiple_accessionIDs = [hdr.split(' ')[0] for hdr in line[1:].split('\x01')]
                                for accessionID in multiple_accessionIDs[:-1]:
                                    out.write(accessionID + '\x01')
                                out.write(multiple_accessionIDs[-1] + '\n')
                            else:
                                accessionID = (line[1:].split(' ')[0])
                                out.write('>' + accessionID + '\n')
                        else:
                            out.write(line)
            logger.info("Fasta database %s has now reduced headers and is saved to %s" % (Path(path_to_db).name, rh_path))

        except FileNotFoundError:
            logger.exception('Database %s not found. No reduction of NCBI database headers.' % path_to_db,
                             exc_info=True)


