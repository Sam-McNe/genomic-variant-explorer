#!/bin/bash

# --- Script Configuration and Robustness ---
set -eEuo pipefail
[[ "${DEBUG:-}" == "1" ]] && set -x

# --- Input and Output ---
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <gene_name>"
    exit 1
fi
gene_name=$1

# Define file paths
GFF_FILE="./Bvulgarisssp_vulgaris_782_EL10.2_2.gene_exons.gff3"
VCF_FILE="./mitchseq_maf1pct_maxmissing95pct_minMeanDP10_maxMeanDP500.vcf.gz"
ANNOTATION_FILE="./Bvulgarisssp_vulgaris_782_EL10.2_2.annotation_info.txt"
OUTPUT_DIR="./output"

# Prepare output directory
mkdir -p "$OUTPUT_DIR"
rm -rf "${OUTPUT_DIR:?}"/*

echo "Processing gene: $gene_name..."

### ======================================================================
### 1. ANNOTATION FILE PROCESSING
### ======================================================================
echo "  -> Extracting annotations..."
ANNOTATION_OUTPUT="${OUTPUT_DIR}/${gene_name}.annotation.txt"

awk -v gene="$gene_name" '
    BEGIN {
        FS=OFS="\t";
        print "Best-hit-arabi-name", "Best-hit-arabi-defline", "Best-hit-clamy-name", "Best-hit-clamy-defline", "Best-hit-rice-name", "Best-hit-rice-defline";
    }
    $1 !~ /^#/ && $2 == gene {
        print $11, $12, $13, $14, $15, $16;
    }
' "$ANNOTATION_FILE" > "$ANNOTATION_OUTPUT"

if [[ $(wc -l < "$ANNOTATION_OUTPUT") -le 1 ]]; then
    echo "Warning: No annotation found for '$gene_name' in column 2 of annotation file." >&2
fi


### ======================================================================
### 2. GFF "MAPPING"
### ======================================================================
echo "  -> Mapping gene coordinates from GFF..."
TEMP_GFF_SLICE="${OUTPUT_DIR}/${gene_name}.temp.gff"
trap 'rm -f "$TEMP_GFF_SLICE"' EXIT

# Create a small GFF slice
awk -v gene="$gene_name" '
    $9 ~ "ID=" gene ";" || 
    $9 ~ "ID=" gene "." ||
    $9 ~ "Name=" gene ";" || 
    $9 ~ "Parent=" gene ";" ||
    $9 ~ "Parent=" gene "."
' "$GFF_FILE" > "$TEMP_GFF_SLICE"

if [[ ! -s "$TEMP_GFF_SLICE" ]]; then
    echo "Error: No GFF records found for '$gene_name'. Check gene name and GFF format (expected ID=<gene>; or Name=<gene>; in column 9)." >&2
    exit 1
fi

# Extract gene-level coordinates
awk '$3 == "gene" {print $1, $4, $5}' "$TEMP_GFF_SLICE" > "${OUTPUT_DIR}/${gene_name}.region.txt"

if [[ ! -s "${OUTPUT_DIR}/${gene_name}.region.txt" ]]; then
    echo "Error: GFF slice for '$gene_name' contains no 'gene' feature in column 3. Features found:" >&2
    awk '{print $3}' "$TEMP_GFF_SLICE" | sort -u >&2
    exit 1
fi

# Extract all features
awk '{print $3,$4,$5}' "$TEMP_GFF_SLICE" > "${OUTPUT_DIR}/${gene_name}.features.txt"

if [[ ! -s "${OUTPUT_DIR}/${gene_name}.features.txt" ]]; then
    echo "Error: Could not extract features from GFF slice for '$gene_name'." >&2
    exit 1
fi

# Read coordinates into variables for the VCF slicing step
read TMPGENECHR TMPSTART TMPEND < "${OUTPUT_DIR}/${gene_name}.region.txt"
if [[ -z "$TMPGENECHR" || -z "$TMPSTART" || -z "$TMPEND" ]]; then
    echo "Error: Could not parse coordinates for '$gene_name'. region.txt contents:" >&2
    cat "${OUTPUT_DIR}/${gene_name}.region.txt" >&2
    exit 1
fi


### ======================================================================
### 3. VCF SLICING
### ======================================================================
echo "  -> Slicing VCF with tabix..."
/project/sugarbeet/samuel/tabix_env/bin/tabix -h "$VCF_FILE" "${TMPGENECHR}:${TMPSTART}-${TMPEND}" > "${OUTPUT_DIR}/${gene_name}.vcf"

if [[ ! -s "${OUTPUT_DIR}/${gene_name}.vcf" ]]; then
    echo "Error: tabix produced an empty VCF for '$gene_name' at ${TMPGENECHR}:${TMPSTART}-${TMPEND}. Check that the chromosome name matches the VCF." >&2
    exit 1
fi

VARIANT_COUNT=$(grep -v "^#" "${OUTPUT_DIR}/${gene_name}.vcf" | wc -l || true)
if [[ "$VARIANT_COUNT" -eq 0 ]]; then
    echo "Warning: VCF for '$gene_name' contains no variants in region ${TMPGENECHR}:${TMPSTART}-${TMPEND}." >&2
fi

echo "Successfully processed $gene_name ($VARIANT_COUNT variants). Files saved to $OUTPUT_DIR"