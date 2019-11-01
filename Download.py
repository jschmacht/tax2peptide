import hashlib
import urllib.request
import logging.config
import wget
from pathlib import Path
from shutil import which

logger = logging.getLogger(__name__)


class Download:

    def __init__(self, url=None, path=None, md5=None):
        """
        :param url: url address for download
        :param path: path to location of partial/full downloaded file
        :param md5: correct md5 hash value
        """
        self.url = url
        self.hash = md5
        self.path = Path(path)
        self.num_bytes = 0
    
	# Validates a file against MD5 hash value
    def validate_file(self):
        """
            return: True if md5 hash value of downloaded file is correct
            """
        calc_md5 = hashlib.md5()
        if self.hash is not None:
            try:
                with open(str(self.path), "rb") as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        calc_md5.update(chunk)
                logger.info("Md5 hash value of correct file: " + str(self.hash))
                logger.info("Md5 hash value of actual file: " + str(calc_md5.hexdigest()))
                return calc_md5.hexdigest() == self.hash
            except FileNotFoundError:
                logger.exception("Path to file to be validated does not exist. No md5 validation possible.", exc_info=True)
        else:
            logging.error("Md5-hash could not be downloaded. No md5 check performed.")
            return False

    # download with python package urllib.request
    def download_with_urllib(self):
        logger.info("Start download from " + self.url + "to location" + str(self.path) + "with urllib.")
        try:
            urllib.request.urlretrieve(self.url, str(self.path))
            logger.info("Download of " + str(self.url) + " file to " + str(self.path) + " was successful.")
        except urllib.error.URLError:
            logger.exception("Download was not successful. Check internet connection and run again. ", exc_info=True)
            exit(1)
        if self.path.stat().st_size < file_size:
            logger.info("Download of the file from " + str(self.path) + " failed. Please run again.")
            exit(1)

    # download with python package urllib.request, resumable if header contains range
    def download_with_range(self, file_size, exist_size):
        """
        :param file-size: size of the complete file
        :param exist-size: size of the partial downloaded file
        """
        logger.info("Start resumable download from " + self.url + "to location" + str(self.path) + ".")
        if self.path.exists():
            # Open file for appending in binary format
            output_file = open(str(self.path), "ab")
            # If the file exists, then download only the remainder
        else:
            output_file = open(str(self.path), "wb")
        header = {"Range": "bytes=%s-" % exist_size}
        req = urllib.request.Request(self.url, headers=header)
        response = urllib.request.urlopen(req)
        n = 0
        while file_size > exist_size and n < 5:
            n += 1
            logger.debug("Continue downloading form url ", self.url, " ", file_size - exist_size, " more bytes")
            while True:
                # read chunks of size 8kb in data and add them to output_file
                data = response.read(8 * 1024)
                if not data:
                    break
                output_file.write(data)
                self.num_bytes += len(data)
        response.close()
        output_file.close()
        logger.debug("downloaded", self.num_bytes, "bytes from", self.url)

    # download with python package wget
    def download_with_wget(self, file_size):
        """
        :param file-size: size of the complete file
        """
        logger.info("Start download from " + self.url + " to location " + str(self.path) + " with wget.")
        try:
            for i in range(5):
                wget.download(str(self.url), str(self.path))
                if self.path.stat().st_size == file_size:
                    logger.info("Download of the file to folder from" + self.url + " to " + str(self.path) + " was successful")
                    break
        # catch error 403: Forbiddden, some server do not accept wget
        except Exception:
            logger.exception("Download of the file from " + self.url + " with wget failed. Switch to urllib download.", exc_info=True)
            self.download_with_urllib()
        if self.path.stat().st_size == 0:
            logger.info("Download of the file from " + self.url + " with wget failed. Switch to urllib download.")
            self.download_with_urllib()
        elif self.path.stat().st_size < file_size:
            logger.info("Download of the file from " + str(self.path) + " is not complete. Please run again.")
            exit(1)

    def download(self):
        try:
            # check to quit the download if the whole file is already downloaded.
            exist_size = self.path.stat().st_size if self.path.exists() else 0
            file_size = int((urllib.request.urlopen(self.url)).headers['Content-Length'])

            if file_size == exist_size:
                logger.info('File ' + str(self.path) + ' was already downloaded from URL ' + self.url + ' .')
            if file_size > exist_size:
                # check if website accept range in header
                if urllib.request.urlopen(self.url).headers['Range']:
                    self.download_with_range(file_size, exist_size)
                # download with wget
                else:
                    self.download_with_wget(file_size)
        except urllib.error.URLError:
            logger.exception('Download from url %s failed.' % self.url, exc_info=True)

        # if md5 hash is given, downloaded file is validated against it
        if self.hash:
            if self.validate_file():
                logger.info('File integrity is md5 checked.')
            else:
                logger.error('Md5 validation of the file %s failed, program continues without md5-check. '
                             % self.url.split('/')[-1])

