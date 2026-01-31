# Genomic Variant Explorer
Interactive Shiny application for querying and visualizing genomic variants across sugar beet germplasm accessions. Developed for USDA Agricultural Research Service to support crop disease resistance research.

## Features

Interactive variant search - Query variants by gene name, genomic position, and accession
Multi-year visualization - Explore variant distributions across multiple years and field trials
Optimized data processing - Backend bash scripts efficiently process large genomic files on-the-fly
Relational data integration - Connects VCF, GFF3, and phenotypic data through common keys
Export functionality - Download filtered results for downstream analysis

## Technical Stack

- R/Shiny - Interactive web application framework
- Bash scripting - Optimized backend data processing
- Docker - Containerized for reproducibility and deployment
- Git - Version control

## Architecture
The application follows a three-tier structure:

1. Data layer: VCF/vcf.GZ/vcf.gz.tbi (variants), GFF3 (gene annotations), and annotation file for best homologs
2. Processing layer: Bash scripts perform on-demand joins and filtering based on user queries
3. Presentation layer: Shiny interface for dynamic querying and visualization

When users submit a search query, the backend scripts:

- Filter variants by genomic coordinates and gene names
- Join across multiple data sources using shared accession IDs
- Return results to the Shiny app for interactive visualization

## Note on Data
This tool was developed for USDA Agricultural Research Service germplasm research. Example data is not included in this repository as it represents ongoing research datasets. The repository demonstrates the application architecture, code structure, and data processing approach.
To use this tool with your own data, you would need:

- VCF file with genomic variants
- GFF3 file with gene annotations
- Phenotype/metadata file with accession information
- Modfiy the existing scripts and app to fit your data

## Development Context
Built as part of genomic data analysis workflows at USDA ARS Fort Collins, focusing on crop disease resistance research and supporting international collaborative breeding programs.

## Author
Sam McNeill\
Biological Data Scientist - Contractor, USDA Agricultural Research Service\
MS Applied Statistics/Data Science, Colorado State University
