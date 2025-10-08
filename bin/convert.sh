#!/usr/bin/bash
# set -euo pipefail
SCRIPT_DIR=$(cd $(dirname $0); pwd)
CONFIG_DIR=$SCRIPT_DIR/../config

# 分割処理の最大行数
SPLIT_THRESHOLD=20000000

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <dir1> [dir2 ...]" >&2
    exit 1
fi

# Turtle ファイルを分割して rapper で処理する関数
process_turtle_file() {
    local file=$1
    echo "Processing $file..." >&2

    # ファイルの行数を取得
    local line_count=$(wc -l < "$file")
    echo "$file contains $line_count lines" >&2

    if [ "$line_count" -gt "$SPLIT_THRESHOLD" ]; then
        echo "File is large, splitting into chunks of $SPLIT_THRESHOLD lines" >&2

        # 一時ディレクトリを作成
        local tmp_dir="tmp_split_$(basename "$file" .ttl)"
        mkdir -p "$tmp_dir"

        # ファイルを分割（トリプルの境界で）
        awk -v threshold="$SPLIT_THRESHOLD" -v outdir="$tmp_dir" '
        BEGIN {
            chunk_num = 0
            line_count = 0
            current_file = sprintf("%s/chunk_%02d", outdir, chunk_num)
        }
        {
            print $0 > current_file
            line_count++
            # `.` で終わる行かつ閾値を超えた場合に次のファイルへ
            if ($0 ~ /\.$/ && line_count >= threshold) {
                close(current_file)
                chunk_num++
                current_file = sprintf("%s/chunk_%02d", outdir, chunk_num)
                line_count = 0
            }
        }
        END {
            close(current_file)
        }
        ' "$file"


        # 分割されたファイルを処理して結合
        > "${file}.processed"

        # プレフィックス部分を保存（通常ファイルの先頭部分）
        grep -E "^@prefix|^@base" "$file" > "${tmp_dir}/prefixes.ttl"

        for chunk in "$tmp_dir"/chunk_*; do
            echo "Processing chunk $chunk" >&2

            # プレフィックスをチャンクの先頭に追加して rapper で処理
            cat "${tmp_dir}/prefixes.ttl" "$chunk" > "${chunk}.with_prefix"
            rapper -i turtle -o turtle "${chunk}.with_prefix" > "${chunk}.processed"

            # プレフィックス部分を除去して結合（最初のチャンクを除く）
            if [ "$chunk" = "$tmp_dir/chunk_aa" ]; then
                cat "${chunk}.processed" >> "${file}.processed"
            else
                grep -v -E "^@prefix|^@base" "${chunk}.processed" >> "${file}.processed"
            fi
        done

        # 元のファイルを置き換え
        mv "${file}.processed" "$file"

        # 一時ディレクトリを削除
        rm -rf "$tmp_dir"
    else
        # サイズが閾値以下なら通常処理
        rapper -i turtle -o turtle "$file" > "${file}.rapper.ttl"
        mv "${file}.rapper.ttl" "$file"
    fi
}

for d in "$@"; do
    (
        set -euo pipefail
        if [ ! -d "$d" ]; then
            echo "Warning: '$d' is not a directory or does not exist, skipping." >&2
            exit 1
        fi

        echo "$d" >&2
        cd "$d"
        python3 $SCRIPT_DIR/rdf_converter_ensembl_db.py $CONFIG_DIR/dbinfo.json

        for f in gene transcript translation exon exon_transcript xref; do
            if [ -f "$f.ttl" ]; then
                process_turtle_file "$f.ttl"
                gzip -f "$f.ttl"
            else
                echo "Warning: $f.ttl not found" >&2
            fi
        done
    )
done
