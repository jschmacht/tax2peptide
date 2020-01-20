import logging
from pathlib import Path
from datetime import date

logger = logging.getLogger(__name__)

class Output:
    """
    creates, modifies, and deletes files & folder, output stores tuples (rel name, abs name, dir|file)
    """
    # question to user, answer y (yes) or n (no)
    @staticmethod
    def ask_user(question, yes="y", no="n", sense_case=False):
        """
        :param question: question to prompt
        :param yes: confirmation input
        :param no: decline input
        :param sense_case: case-sensitivity (True/False)
        :return: user prompt returns True (user confirmed) or False (user did not confirm)
        """
        while True:
            prompt = input(question + " [" + str(yes) + "/" + str(no) + "] ")
            if not sense_case:
                prompt = prompt.lower()
                yes = yes.lower()
                no = no.lower()
            if prompt == yes:
                return True
            elif prompt == no:
                return False

    # check is file exist and ask user if they want to override this file
    @staticmethod
    def check_file(path_to_file):
        if path_to_file.is_file():
            if Output.ask_user("Do you want to delete %s and download it again?" % path_to_file.name):
                try:
                    path_to_file.unlink()
                    logger.info("File %s deleted." % path_to_file)
                    return False
                except OSError:
                    logger.exception('Can not remove file %s.'% path_to_file, exc_info=True)
            else:
                return True
        return False

    # check which files exist in database folder
    @staticmethod
    def check_content(file_list):
        """
        :param file_list: contains directions of files to be checked
        :return booleans
        """
        return_list = []
        for file in file_list:
            return_list.append(Output.check_file(file))
        return tuple(return_list)

    @staticmethod
    def check_files_exist(file_list):
        """
                :param file_list: contains directions of files to be checked
                :return list of bools
        """
        return_list = []
        for file in file_list:
            file = Path(file)
            if file.is_file():
                return_list.append(True)
            else:
                return_list.append(False)
        return tuple(return_list)

    # create new direction for the taxon_specific database (taxon_specified_db_DATE) , default name= taxon_database.fasta
    @staticmethod
    def createDir(path_cwd):
        """
        :param path_cwd: new path starts in working direction
                """
        output_path = path_cwd / ('taxon_specified_db' + '_' + str(date.today())) / 'taxon_specific_database.fasta'
        if not output_path.parents[0].exists():
            try:
                output_path.parents[0].mkdir()
            except OSError:
                print("Error creating new direction for saving new generated taxon specified peptide database.")
                exit(1)
        number = 1
        if output_path.is_file() or (output_path.parents[0] / (output_path.stem + '_nr.fasta')).is_file():
            output_path = Path.cwd() / ('taxon_specified_db' + '_' + str(date.today())) / \
                          ('taxon_specific_database' + '_' + str(number) + '.fasta')
            while output_path.is_file() or (output_path.parents[0] / (output_path.stem + '_nr.fasta')).is_file():
                number += 1
                output_path = Path.cwd() / ('taxon_specified_db' + '_' + str(date.today())) / \
                              ('taxon_specific_database' + '_' + str(number) + '.fasta')
        return output_path
