#!/usr/bin/env bash
# Step 2: download the 53 paired-end samples listed in download_manifest.tsv
# from ENA (~3.3 GB) and give them the Casava-style file names QIIME2 requires:
#
#     <sample-id>_S<n>_L001_R[12]_001.fastq.gz
#
# Every file is checked against the MD5 ENA publishes. Safe to re-run: files
# already present with the right MD5 are skipped, anything else is
# re-downloaded.
#
# Usage: 02_download_reads.sh [reads_dir]        (default: ./raw_reads)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="$HERE/download_manifest.tsv"
READS="${1:-raw_reads}"
mkdir -p "$READS"

fetch () {  # fetch <url> <md5> <destination>
    local url="$1" md5="$2" dest="$3"
    if [ -s "$dest" ] && [ "$(md5sum < "$dest" | cut -d' ' -f1)" = "$md5" ]; then
        echo "    ok (already downloaded): $(basename "$dest")"
        return
    fi
    curl --fail --silent --show-error --location --retry 5 --retry-delay 10 -o "$dest.part" "https://$url"
    if [ "$(md5sum < "$dest.part" | cut -d' ' -f1)" != "$md5" ]; then
        rm -f "$dest.part"
        echo "MD5 mismatch for $url" >&2
        exit 1
    fi
    mv "$dest.part" "$dest"
    echo "    ok: $(basename "$dest")"
}

n=0
while IFS=$'\t' read -r sample_id run_accession fastq_ftp fastq_md5; do
    n=$((n + 1))
    echo "[$n] $sample_id ($run_accession)"
    fetch "${fastq_ftp%%;*}" "${fastq_md5%%;*}" "$READS/${sample_id}_S${n}_L001_R1_001.fastq.gz"
    fetch "${fastq_ftp##*;}" "${fastq_md5##*;}" "$READS/${sample_id}_S${n}_L001_R2_001.fastq.gz"
done < <(tail -n +2 "$MANIFEST")

echo "$n samples in $READS ($(du -sh "$READS" | cut -f1))"
