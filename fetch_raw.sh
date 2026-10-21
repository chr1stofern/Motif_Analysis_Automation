#!/usr/bin/env bash
#
# fetch_raw.sh — download R1 FASTQs for the deposited Cas13 dataset.
#
# Reads samples.csv (comma-separated; cols: accession, replicate, protein, type, ...),
# writes {raw_dir}/{replicate}_{protein}-{type}.fastq.gz (R1 only).
#
# - Rows with an empty 'accession' are skipped (local-data users).
# - For each accession, ENA's API is queried for the real R1 path.
# - Files already present (non-empty) are skipped (idempotent).
# - If a run is not mirrored to ENA, sracha (if available) fetches it from NCBI.
# - R2 (mate) is never kept.
#
# Requires: curl  (and sracha for runs not on ENA)
#
# Usage: ./fetch_raw.sh [samples.csv] [raw_dir]
#
set -euo pipefail

SAMPLES_CSV="${1:-samples.csv}"
RAW_DIR="${2:-raw}"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

[[ -f "$SAMPLES_CSV" ]] || { echo "ERROR: $SAMPLES_CSV not found" >&2; exit 2; }
command -v curl >/dev/null || { echo "ERROR: curl not on PATH" >&2; exit 2; }

mkdir -p "$RAW_DIR"
# Ensure the CSV ends with a newline so the last sample is never skipped
[[ $(tail -c1 "$SAMPLES_CSV" | wc -l) -eq 0 ]] && echo >> "$SAMPLES_CSV"

# Clean up leftover partial / temporary files from previous interrupted runs
echo "Cleaning leftover partial files in $RAW_DIR ..."
find "$RAW_DIR" -maxdepth 1 -type f \( \
    -name '*.part' \
    -o -name 'SRR*_1.fastq.gz' \
    -o -name 'SRR*_2.fastq.gz' \
\) -print -delete 2>/dev/null || true

tail -n +2 "$SAMPLES_CSV" | tr -d '\r' | while IFS=, read -r srr rep prot type _; do
    [[ -z "$srr" ]] && { echo "skip ${prot}-${type} (no accession)"; continue; }

    out="$RAW_DIR/${rep}_${prot}-${type}.fastq.gz"
    [[ -s "$out" ]] && { echo "skip $out (exists)"; continue; }

    # Retry ENA a few times (handles transient empty replies)
    ftp_path=""
    for try in 1 2 3; do
        ftp_path=$(curl -s --connect-timeout 15 \
            "${ENA_API}?accession=${srr}&result=read_run&fields=fastq_ftp" \
            | tail -1 | cut -f2 | cut -d';' -f1)
        [[ -n "$ftp_path" ]] && break
        sleep $((try * 2))
    done

    if [[ -n "$ftp_path" ]]; then
        echo "fetch $srr (ENA) -> $out"
        if curl -sSL --fail --retry 5 --retry-delay 5 "https://${ftp_path}" -o "$out.part"; then
            mv "$out.part" "$out"
        else
            echo "FAIL $srr (ENA download)" >&2
            rm -f "$out.part"
            continue
        fi
    else
        echo "MISSING $srr (${prot}-${type}) — trying sracha fallback"
        if ! command -v sracha >/dev/null; then
            echo "FAIL $srr (not on ENA; sracha not on PATH)" >&2
            continue
        fi
        if ( cd "$RAW_DIR" && sracha get "$srr" ); then
            mv "$RAW_DIR/${srr}_1.fastq.gz" "$out"
            rm -f "$RAW_DIR/${srr}_2.fastq.gz"
            echo "fetch $srr (sracha) -> $out"
        else
            echo "FAIL $srr (sracha)" >&2
            rm -f "$RAW_DIR/${srr}_1.fastq.gz" "$RAW_DIR/${srr}_2.fastq.gz"
            continue
        fi
    fi

    # Quick integrity check
    if ! gzip -t "$out" 2>/dev/null; then
        echo "CORRUPT $out – removing so it can be retried" >&2
        rm -f "$out"
    fi
done
echo
echo "----- summary -----"
present=$(find "$RAW_DIR" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)
echo "FASTQs in $RAW_DIR: $present"
echo "For the deposited dataset, expect 22 files."
echo "If fewer, check the README's 'Handling missing runs' section."