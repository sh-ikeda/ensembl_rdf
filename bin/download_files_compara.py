import os
import sys
import json
import datetime
from ftplib import FTP


def download_files(ftp, directory, dbs):
    start_dir = ftp.pwd()
    ftp.cwd(directory)
    files = ftp.nlst()
    os.makedirs(directory, exist_ok=True)
    for file in files:
        if file in dbs:
            dt_now = datetime.datetime.now()
            print(f"[{dt_now}] Downloading: {directory}/{file}", file=sys.stderr)
            with open(directory + "/" + file, "wb") as f:
                ftp.retrbinary("RETR %s" % file, f.write)
    ftp.cwd(start_dir)


def process_directory(ftp, directories, dbs):
    for directory in directories:
        ftp.cwd(directory)

        subdirectories = ftp.nlst()
        for subdirectory in subdirectories:
            if "_compara_" in subdirectory:
                dt_now = datetime.datetime.now()
                print(f"[{dt_now}] {subdirectory}", file=sys.stderr)
                download_files(ftp, subdirectory, dbs)


def main():
    input_dbinfo_file = sys.argv[1]
    with open(input_dbinfo_file, "r") as f:
        dbinfo = json.load(f)
    dbs = [dbinfo[k]["filename"] for k in dbinfo]

    ## vertebrate
    ftp_url = "ftp.ensembl.org"
    ftp_dirs = ["/pub/current_mysql/"]

    ftp = FTP(ftp_url)
    ftp.login()
    process_directory(ftp, ftp_dirs, dbs)
    ftp.quit()

    ## genomes
    ftp_url = "ftp.ensemblgenomes.ebi.ac.uk"
    genomes_sections = ["metazoa", "plants", "protists", "fungi"]
    ftp_dirs = ["/pub/current/"+s+"/mysql/" for s in genomes_sections]
    ftp_dirs.append("/pub/pan_ensembl/current/mysql/")

    ftp = FTP(ftp_url)
    ftp.login()
    process_directory(ftp, ftp_dirs, dbs)
    ftp.quit()

if __name__ == "__main__":
    main()
