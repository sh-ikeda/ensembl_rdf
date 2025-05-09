# Ensembl RDF converter by DBCLS
## Data download
```
$ python3 /path/to/ensembl_rdf/bin/download_files.py ftp.ensembl.org /pub/current_mysql/
```
## RDF conversion
In the directory where `download_files.py` were executed, run `convert.sh`.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh
```
Options:  
- `-a`: Execute only for the `acanthochromis_polyacanthus_core_` dataset. Use for a minimal test.
- `-t`: Execute only for the `homo_sapiens_core_` dataset.

## RDF schema
[RDF-config](https://github.com/dbcls/rdf-config/blob/master/config/ensembl/model.yaml)
