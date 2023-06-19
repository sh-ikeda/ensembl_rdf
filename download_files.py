import os
import sys
import json
import datetime
from ftplib import FTP


def download_files(ftp, directory):
    start_dir = ftp.pwd()
    ftp.cwd(directory)
    files = ftp.nlst()
    # print(files)
    os.makedirs(directory, exist_ok=True)
    for file in files:
        if file in dbs:
            dt_now = datetime.datetime.now()
            print(f'[{dt_now}] Downloading: {directory}/{file}', file=sys.stderr)
            with open(directory + "/" + file, 'wb') as f:
                ftp.retrbinary('RETR %s' % file, f.write)
    ftp.cwd(start_dir)


def process_directory(ftp, directory):
    ftp.cwd(directory)

    subdirectories = ftp.nlst()
    for subdirectory in subdirectories:
        if 'musculus_core_' in subdirectory:
            dt_now = datetime.datetime.now()
            print(f'[{dt_now}] {subdirectory}', file=sys.stderr)
            download_files(ftp, subdirectory)


BASE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"
with open(BASE_DIR+"dbinfo.json", "r") as f:
    dbinfo = json.load(f)
dbs = [dbinfo[k]["filename"] for k in dbinfo]
print(dbs)
ftp = FTP('ftp.ensembl.org')
ftp.login()
process_directory(ftp, '/pub/current_mysql/')
ftp.quit()
