
import gzip
import shutil
import requests
import os
from constants import DATA_DIRECTORY, STRING_URL, LIT_BM, HURI_URL, HI_UNION, HUMAN_TF, STRING_ALIAS, \
        STRING_URL_PATH, STRING_ALIAS_PATH, HURI_URL_PATH, HI_UNION_PATH, LIT_BM_PATH, HUMAN_TF_PATH, \
        STRING_PROTEIN_LIST_URL, STRING_PROTEIN_LIST_URL_PATH, STRING_BASE_PATH, STRING_EXTRACTED_ALIAS_PATH, \
        STRING_EXTRACTED_URL_PATH, STRING_EXTRACTED_PROTEIN_LIST_URL_PATH, HURI_BASE_PATH, HUMAN_TF_BASE_PATH, HWG_BASE_PATH, \
        HUMAN_TF_IDLIST_URL, HUMAN_TF_IDLIST_PATH

class Downloader():
    '''
    Basic downloader class for downloading files from url to a directory
    '''
    def __init__(self, url: str, destination: str, extractpath: str = None):
        self.url = url
        self.destination = destination
        self.extractpath = extractpath

    def download(self):
        '''
        Function to download the file from url as chunks and save it to destination.
        '''
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)

        print(f'Starting download from {self.url} to {self.destination}')
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)

        filename = self.url.split('/')[-1].replace(" ", "_")
        filepath = os.path.join(self.destination, filename)

        r = requests.get(self.url, stream=True)
        if r.status_code == 200:
            with open(filepath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rb') as f_in:
                with open(self.extractpath, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"Downloaded {filename} to {self.destination}. It is a gzipped file.")


def setup_db():

    downloader = Downloader(url=STRING_ALIAS, destination=STRING_BASE_PATH, extractpath=STRING_EXTRACTED_ALIAS_PATH)
    downloader.download()
    downloader = Downloader(url=STRING_URL, destination=STRING_BASE_PATH, extractpath = STRING_EXTRACTED_URL_PATH)
    downloader.download()
    downloader = Downloader(url=STRING_PROTEIN_LIST_URL, destination=STRING_BASE_PATH, extractpath=STRING_EXTRACTED_PROTEIN_LIST_URL_PATH)
    downloader.download()
    downloader = Downloader(url=HURI_URL, destination=HURI_BASE_PATH)
    downloader.download()
    downloader = Downloader(url=HI_UNION, destination=HURI_BASE_PATH)
    downloader.download()
    downloader = Downloader(url=LIT_BM, destination=HURI_BASE_PATH)
    downloader.download()
    downloader = Downloader(url=HUMAN_TF, destination=HUMAN_TF_BASE_PATH)
    downloader.download()
    downloader = Downloader(url=HUMAN_TF_IDLIST_URL, destination=HUMAN_TF_BASE_PATH)
    downloader.download()
    



if __name__ == "__main__":
    setup_db()
