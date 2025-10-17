# Ensembl RDF converter by DBCLS
## Data download
```
$ python3 /path/to/ensembl_rdf/bin/download_files.py ftp.ensembl.org /pub/current_mysql/
```
## RDF conversion
In the directory where `download_files.py` were executed, run `convert.sh`.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh *_core_*
```

If only some types of output are needed, the `-t` option can be specified.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh *_core_* -t gene -t transcript
```

Note that the direct use of the python script needs only one `-t` to specify multiple types.
```
$ python3 rdf_converter_ensembl_db.py config/dbinfo.json path/to/input/ -t gene transcript
```

## RDF schema
[RDF-config](https://github.com/dbcls/rdf-config/blob/master/config/ensembl/model.yaml)
