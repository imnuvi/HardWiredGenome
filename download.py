
import gzip
import shutil
import requests
import zipfile

import os
from constants import DATA_DIRECTORY, STRING_URL, LIT_BM, HURI_URL, HI_UNION, HUMAN_TF, STRING_ALIAS, \
        STRING_URL_PATH, STRING_ALIAS_PATH, HURI_URL_PATH, HI_UNION_PATH, LIT_BM_PATH, HUMAN_TF_PATH, \
        STRING_PROTEIN_LIST_URL, STRING_PROTEIN_LIST_URL_PATH, STRING_BASE_PATH, STRING_EXTRACTED_ALIAS_PATH, \
        STRING_EXTRACTED_URL_PATH, STRING_EXTRACTED_PROTEIN_LIST_URL_PATH, HURI_BASE_PATH, HUMAN_TF_BASE_PATH, HWG_BASE_PATH, \
        HUMAN_TF_IDLIST_URL, HUMAN_TF_IDLIST_PATH, HTF_MOTIFS_URL, HTF_MOTIFS_DIR, HTF_MOTIFS_PATH, REFERENCE_GTF_HG38_URL, REFERENCE_GTF_HG38_PATH, REFERENCE_DIR, \
        REFERENCE_GENES_BED_URL, REFERENCE_GENES_BED_PATH, REFERENCE_CHROM_SIZES_PATH, REFERENCE_CHROM_SIZES_URL, REFERENCE_GENOME_URL, REFERENCE_GENOME_PATH, \
        GENCODE_ANNOTATION_URL, GENCODE_ANNOTATION_PATH, REFERENCE_FNA_URL, REFERENCE_FNA_PATH, REFERENCE_GENOME_GENCODE_URL, REFERENCE_GENOME_GENCODE_PATH

class Downloader():
    '''
    Basic downloader class for downloading files from url to a directory
    '''
    def __init__(self, url: str, destination: str, extractpath: str = None, filepath: str = None):
        self.url = url
        self.destination = destination
        self.extractpath = extractpath
        self.filepath = filepath

    def download(self):
        '''
        Function to download the file from url as chunks and save it to destination.
        '''
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)

        print(f'Starting download from {self.url} to {self.destination}')
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)


        if self.filepath:
            filepath = self.filepath
            filename = filepath.split('/')[-1]
        else:
            filename = self.url.split('/')[-1].replace(" ", "_")
            filepath = os.path.join(self.destination, filename)

        r = requests.get(self.url, stream=True)
        if r.status_code == 200:
            with open(filepath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

        print('File downloaded and written to disk')
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rb') as f_in:
                with open(self.extractpath, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"Downloaded {filename} to {self.destination}. It is a gzipped file.")

        if filepath.endswith('.zip'):
            with zipfile.ZipFile(filepath, 'r') as zip_ref:
                zip_ref.extractall(self.extractpath)
            print(f"Downloaded {filename} to {self.destination}. It is a zipped file.")

        print(f"Finished Downloading {filename} to {filepath}")
    def queried_file(self):
        pass


def setup_db():

    # downloader = Downloader(url=STRING_ALIAS, destination=STRING_BASE_PATH, extractpath=STRING_EXTRACTED_ALIAS_PATH)
    # downloader.download()
    # downloader = Downloader(url=STRING_URL, destination=STRING_BASE_PATH, extractpath = STRING_EXTRACTED_URL_PATH)
    # downloader.download()
    # downloader = Downloader(url=STRING_PROTEIN_LIST_URL, destination=STRING_BASE_PATH, extractpath=STRING_EXTRACTED_PROTEIN_LIST_URL_PATH)
    # downloader.download()
    # downloader = Downloader(url=HURI_URL, destination=HURI_BASE_PATH)
    # downloader.download()
    # downloader = Downloader(url=HI_UNION, destination=HURI_BASE_PATH)
    # downloader.download()
    # downloader = Downloader(url=LIT_BM, destination=HURI_BASE_PATH)
    # downloader.download()
    # downloader = Downloader(url=HUMAN_TF, destination=HUMAN_TF_BASE_PATH)
    # downloader.download()
    # downloader = Downloader(url=HUMAN_TF_IDLIST_URL, destination=HUMAN_TF_BASE_PATH)
    # downloader.download()
    # downloader = Downloader(url=HTF_MOTIFS_URL, destination=HTF_MOTIFS_PATH, extractpath=HTF_MOTIFS_DIR)
    # downloader.download()  
    # downloader = Downloader(url=REFERENCE_GTF_HG38_URL, destination=REFERENCE_DIR)
    # downloader.download()
    # downloader = Downloader(url=REFERENCE_GENES_BED_URL, destination=REFERENCE_DIR, filepath=REFERENCE_GENES_BED_PATH)
    # downloader.download()
    # downloader = Downloader(url=REFERENCE_CHROM_SIZES_URL, destination=REFERENCE_DIR)
    # downloader.download()
    # downloader = Downloader(url=REFERENCE_GENOME_URL, destination=REFERENCE_DIR, extractpath=REFERENCE_GENOME_PATH)
    # downloader.download()
    # downloader = Downloader(url=GENCODE_ANNOTATION_URL, destination=REFERENCE_DIR, extractpath=GENCODE_ANNOTATION_PATH)
    # downloader.download()
    # downloader = Downloader(url=REFERENCE_FNA_URL, destination=REFERENCE_DIR, extractpath=REFERENCE_FNA_PATH)
    # downloader.download()
    downloader = Downloader(url=REFERENCE_GENOME_GENCODE_URL, destination=REFERENCE_DIR, extractpath=REFERENCE_GENOME_GENCODE_PATH)
    downloader.download()
    
    



if __name__ == "__main__":
    setup_db()
