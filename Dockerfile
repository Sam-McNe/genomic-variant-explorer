# Start from a base image that has R and Shiny pre-installed
FROM rocker/shiny:4.4.1

# Install system tools (tabix, bcftools, etc.)
RUN apt-get update && apt-get install -y \    
	tabix \
	&& rm -rf /var/lib/apt/lists/*

# Install R packages your app needs
RUN R -e "install.packages(c('shiny', 'ggplot2', 'dplyr', 'DT', 'vcfR'), repos='https://cloud.r-project.org/')"

# Create app directory, data directory (for input files), and output directory (for generated files)
RUN mkdir -p /app/data /app/output

# Copy your files into the container
COPY app.R /app/
COPY variant_mining_tool_1.sh /app/data
RUN chmod +x /app/data/variant_mining_tool_1.sh

# Copy your data files (VCF, GFF, etc.) to /app/data
COPY *.vcf.gz /app/data
COPY *.annotation_info.txt /app/data
COPY *.gff3 /app/data
COPY *.vcf.gz.tbi /app/data

# Set working directory
WORKDIR /app

# Expose port 3838 for Shiny
EXPOSE 3838

# Run the app
CMD ["R", "-e", "shiny::runApp('/app/app.R', host='0.0.0.0', port=3838)"]
