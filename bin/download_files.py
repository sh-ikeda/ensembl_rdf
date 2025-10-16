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
            print(f"[{dt_now}] Downloading: {directory}/{file}", file=sys.stderr)
            with open(directory + "/" + file, "wb") as f:
                ftp.retrbinary("RETR %s" % file, f.write)
    ftp.cwd(start_dir)


def process_directory(ftp, directory):
    ftp.cwd(directory)

    subdirectories = ftp.nlst()
    for subdirectory in subdirectories:
        if "_core_" in subdirectory:
            dt_now = datetime.datetime.now()
            print(f"[{dt_now}] {subdirectory}", file=sys.stderr)
            download_files(ftp, subdirectory)


BASE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"
CONFIG_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../config/"
ftp_url = sys.argv[1]  # e.g. "ftp.ensembl.org"
ftp_dir = sys.argv[2]  # e.g. "/pub/current_mysql/"

with open(CONFIG_DIR+"dbinfo.json", "r") as f:
    dbinfo = json.load(f)
dbs = [dbinfo[k]["filename"] for k in dbinfo]
ftp = FTP(ftp_url)
ftp.login()
process_directory(ftp, ftp_dir)
ftp.quit()
