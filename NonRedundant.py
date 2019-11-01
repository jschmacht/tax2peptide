import logging.config
from Bio import SeqIO
from pathlib import Path

logger = logging.getLogger(__name__)


class NonRedundant:
    # modified https://biopython.org/wiki/Sequence_Cleaner
    @staticmethod
    def sequence_cleaner(fasta_file, nr_file, with_taxon_ID):
        """
        :param fasta_file: fasta file to make non redundant
        :param nr_file: output file, non redundant
        :param with_taxon_ID: True for uniprot, False for NCBI, different delimiter in header
        """
        sequences = {}
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
            #remove redundant database:
            fasta_file.unlink()
        except OSError:
            logger.exception('Not able to write non redundant database.', exc_info=True)

        logger.info("Fasta database %s is now non redundant " % fasta_file)

