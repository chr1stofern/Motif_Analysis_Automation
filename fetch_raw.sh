#!/usr/bin/env bash
#
# fetch_raw.sh — download R1 FASTQs from ENA for the deposited Cas13 dataset.
# Reads samples.csv (comma-separated; cols: accession, replicate, protein, type, ...),
# writes {raw_dir}/{replicate}_{protein}-{type}.fastq.gz (R1 only).
# Rows with an empty 'accession' are skipped (local-data users).
#
# Usage: ./fetch_raw.sh [samples.csv] [raw_dir]
#
set -euo pipefail

SAMPLES_CSV="${1:-samples.csv}"
RAW_DIR="${2:-raw}"
ENA="https://ftp.sra.ebi.ac.uk/vol1/fastq"

tail -n +2 "$SAMPLES_CSV" | tr -d '\r' | while IFS=, read -r srr rep prot type _; do
    [[ -z "$srr" ]] && { echo "skip ${prot}-${type} (no accession)"; continue; }

    out="$RAW_DIR/${rep}_${prot}-${type}.fastq.gz"
    [[ -s "$out" ]] && { echo "skip $out (exists)"; continue; }

    n="${srr#SRR}"
    url="$ENA/SRR${n:0:3}/$(printf '%03d' $((10#$n % 100)))/$srr/${srr}_1.fastq.gz"
    
    mkdir -p "$RAW_DIR"
    rm -f "$out"
    echo "fetch $srr -> $out"
    if ! curl -sSL --fail "$url" -o "$out"; then
        echo "FAIL $srr" >&2
        rm -f "$out"
    fi
done