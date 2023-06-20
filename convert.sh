#!/usr/bin/bash
set -euo pipefail
SCRIPT_DIR=$(cd $(dirname $0); pwd)

for d in $( ls . | grep "_core_" ); do
    if [ ! -d $d ]; then
        continue
    fi
    cd $d
    echo $d 1>&2
    python3 $SCRIPT_DIR/rdf_converter_ensembl_db.py $SCRIPT_DIR/dbinfo.json

    #echo "Validating turtle files..."
    for f in gene transcript translation exon exon_transcript xref ; do
        rapper -i turtle -o turtle $f.ttl > $f.rapper.ttl
        mv $f.rapper.ttl $f.ttl
        gzip -f $f.ttl
    done
    cd ..
done
