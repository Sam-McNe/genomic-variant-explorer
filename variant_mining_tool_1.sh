#!/bin/bash

gene_name=$1

# Input files are in /app/data
GFF_FILE=$(ls /app/data/*.gff3 | head -n 1)
VCF_FILE=$(ls /app/data/*.vcf.gz | head -n 1)
ANNOTATION_FILE=$(ls /app/data/*.annotation_info.txt | head -n 1)

# Output files go to /app/output
OUTPUT_DIR=/app/output
ANNOTATION_OUTPUT=$OUTPUT_DIR/$gene_name.annotation.txt

###### ANNOTATION FILE ######
# Get info from annotation file
grep -w "$gene_name" "$ANNOTATION_FILE" | awk -v FS='\t' -v OFS='\t' '{print $11, $12, $13, $14, $15, $16}' > $ANNOTATION_OUTPUT

# Create a temporary file with the headers
HEADERS="Best-hit-arabi-name\tBest-hit-arabi-defline\tBest-hit-clamy-name\tBest-hit-clamy-defline\tBest-hit-rice-name\tBest-hit-rice-defline"
echo -e "$HEADERS" > $OUTPUT_DIR/temp_file.txt

# Append and replace temp file
cat $ANNOTATION_OUTPUT >> $OUTPUT_DIR/temp_file.txt
mv $OUTPUT_DIR/temp_file.txt $ANNOTATION_OUTPUT

###### GFF ######

grep -w "$gene_name" "$GFF_FILE" | grep -w "gene" | awk '{print $1, $4, $5}' > $OUTPUT_DIR/$gene_name.region.txt
grep -w "$gene_name" "$GFF_FILE" | awk '{print $3,$4,$5}' > $OUTPUT_DIR/$gene_name.features.txt

read TMPGENECHR TMPSTART TMPEND < $OUTPUT_DIR/$gene_name.region.txt

if [[ -z "$TMPGENECHR" || -z "$TMPSTART" || -z "$TMPEND" ]]; then
    echo "Error: Could not extract gene coordinates for $gene_name"
    exit 1
fi

###### VCF ######

if [[ ! -f "$VCF_FILE" ]] || [[ ! -f "$VCF_FILE.tbi" ]]; then
    echo "Error: VCF file not found or not indexed"
    exit 1
fi

tabix -h "$VCF_FILE" ${TMPGENECHR}:${TMPSTART}-${TMPEND} > $OUTPUT_DIR/$gene_name.vcf

# Add headers to region and features files
# echo -e "chr\tstart\tend" | cat - $OUTPUT_DIR/$gene_name.region.txt > $OUTPUT_DIR/temp && mv $OUTPUT_DIR/temp $OUTPUT_DIR/$gene_name.region.txt
# echo -e "feature_type\tstart\tend" | cat - $OUTPUT_DIR/$gene_name.features.txt > $OUTPUT_DIR/temp && mv $OUTPUT_DIR/temp $OUTPUT_DIR/$gene_name.features.txt

# echo "Successfully processed $gene_name - files saved to $OUTPUT_DIR"
