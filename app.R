# Variant Mining App by Sam McNeill | USDA-ARS | January 2026

# Displays gene information, structure, and variant genotypes with filtering and export


# ===== Download / Load Libraries =====

packages <- c("shiny", "ggplot2", "dplyr", "DT", "vcfR")

missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) install.packages(missing_packages)

invisible(lapply(packages, library, character.only = TRUE))


# ===== Directory Config =====

DATA_DIR <- "/app/"


# ===== Gene Diagram Function =====

plot_gene_structure <- function(features_df, variants_df = NULL, gene_name = "") {
  
  feature_hierarchy <- data.frame(
    feature_type = c("gene", "mRNA", "CDS", "five_prime_UTR", "three_prime_UTR",
                     "exon", "intron"),
    track = c(1, 2, 3, 4, 5, 6, 7),
    color = c("gray60", "gray40", "darkblue", "lightgreen", "lightgreen",
              "steelblue", "gray80"),
    stringsAsFactors = FALSE
  )
  
  plot_data <- features_df %>%
    left_join(feature_hierarchy, by = "feature_type") %>%
    mutate(
      track = ifelse(is.na(track), max(feature_hierarchy$track) + 1, track),
      color = ifelse(is.na(color), "gray50", color)
    )
  
  gene_start <- min(plot_data$start)
  gene_end <- max(plot_data$end)
  
  p <- ggplot(plot_data) +
    geom_rect(aes(xmin = start, xmax = end,
                  ymin = track - 0.35, ymax = track + 0.35,
                  fill = feature_type),
              color = "black", size = 0.3) +
    geom_text(aes(x = gene_start - (gene_end - gene_start) * 0.02,
                  y = track, label = feature_type),
              hjust = 1, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = setNames(plot_data$color, plot_data$feature_type)) +
    scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      expand = expansion(mult = c(0.15, 0.05))
    ) +
    labs(
      title = paste("Gene Structure:", gene_name),
      x = "Genomic Position (bp)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 10)
    )
  
  # Add variants as vertical lines if provided
  
  if (!is.null(variants_df) && nrow(variants_df) > 0) {
    p <- p +
      geom_segment(data = variants_df,
                   aes(x = POS, xend = POS, y = 0.5, yend = max(plot_data$track) + 0.5,
                       color = TYPE),
                   alpha = 0.7, linewidth = 1) +
      scale_color_manual(values = c("SNP" = "red", "INDEL" = "orange", "OTHER" = "gray"),
                         name = "Variant Type")
  }
  
  return(p)
}

# ===== UI =====

ui <- fluidPage(
  
  titlePanel("Variant Mining Tool"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      
      textInput("gene_name", "Enter Gene Name:", placeholder = "e.g., Bevul.2G077600"),
      
      actionButton("submit", "Load Gene", class = "btn-primary", 
                   style = "width: 100%; margin-bottom: 20px;"),
      
      hr(),
      
      h5("Variant Filters"),
      
      selectInput("variant_type", "Variant Type:",
                  choices = c("All" = "all", "SNP" = "SNP", "INDEL" = "INDEL"),
                  selected = "all"),
      
      checkboxInput("coding_only", "Coding variants only", value = FALSE),
      
      hr(),
      
      downloadButton("download_genotypes", "Download Genotypes", 
                     style = "width: 100%; margin-bottom: 10px;"),
      
      downloadButton("download_vcf", "Download VCF", 
                     style = "width: 100%; margin-bottom: 10px;")
    ),
    
    mainPanel(
      width = 9,
      
      # Gene location card
      uiOutput("gene_location_card"),
      
      br(),
      
      # Gene diagram
      plotOutput("gene_diagram", height = "350px"),
      
      br(),
      
      # Two column layout for features and homologs
      fluidRow(
        column(6, 
               wellPanel(
                 h4("Gene Features"),
                 DTOutput("features_table")
               )
        ),
        column(6,
               wellPanel(
                 h4("Top Homologs"),
                 DTOutput("homologs_table")
               )
        )
      ),
      
      br(),
      
      # Variant summary
      wellPanel(
        h4("Variant Summary"),
        uiOutput("variant_summary")
      ),
      
      br(),
      
      # Genotypes table
      conditionalPanel(
        condition = "output.gene_loaded && output.has_variants",
        h4("Genotypes"),
        DTOutput("genotypes_table")
      ),
      conditionalPanel(
        condition = "output.gene_loaded && !output.has_variants",
        h4("Genotypes"),
        div(style = "text-align: center; padding: 40px; color: #999;",
            h4("No variants found in this gene", style = "color: #FF9800;"),
            p("Select a different gene.")
        )
      )
    )
    
    
  )
)

# ===== Server =====

server <- function(input, output, session) {
  
  # Reactive data storage
  gene_data <- reactiveValues(
    location = NULL,
    features = NULL,
    homologs = NULL,
    vcf = NULL,
    variants = NULL,
    current_gene = NULL
  )
  
  # Load data when button clicked
  observeEvent(input$submit, {
    req(input$gene_name)
    
    gene_data$location <- NULL
    gene_data$features <- NULL
    gene_data$homologs <- NULL
    gene_data$vcf <- NULL
    gene_data$variants <- NULL
    gene_data$current_gene <- input$gene_name
    
    tryCatch({
      # Show progress notification
      showNotification("Gathering data...", type = "message", duration = NULL, id = "generating")
      
      result <- system2("bash", 
                        args = c("./data/variant_mining_tool_1.sh", gene_data$current_gene),
                        stdout = TRUE,
                        stderr = TRUE,
                        wait = TRUE)
        
      
      removeNotification(id = "generating")
      
      # Check if files were created
      location_file <- file.path(DATA_DIR, "output", paste0(gene_data$current_gene, ".region.txt"))
      if (!file.exists(location_file)) {
        stop("Data files could not be found.")
      }
      
      # Read location file (chr, start, end)
      gene_data$location <- read.table(location_file, header = FALSE, 
                                       col.names = c("chr", "start", "end"))
      
      # Read features file (feature_type, start, end)
      features_file <- file.path(DATA_DIR, "output", paste0(gene_data$current_gene, ".features.txt"))
      gene_data$features <- read.table(features_file, header = FALSE, 
                                       col.names = c("feature_type", "start", "end"))
      
      # Read homologs file (6 columns, 2 rows)
      homologs_file <- file.path(DATA_DIR, "output", paste0(gene_data$current_gene, ".annotation.txt"))
      gene_data$homologs <- read.table(homologs_file, header = TRUE, sep = "\t")
      colnames(gene_data$homologs) <- c("Best hit Arabidopsis name","Best hit Arabidopsis definition",
                                        "Best hit Clamydomonas name","Best hit Clamydomonas definition",
                                        "Best hit Rice name","Best hit Rice definition")
      
      # Read VCF file
      vcf_file <- file.path(DATA_DIR, "output", paste0(gene_data$current_gene, ".vcf"))
      gene_data$vcf <- read.vcfR(vcf_file, verbose = FALSE)
      
      # Extract variant info
      fix_data <- getFIX(gene_data$vcf)
      
      # Handling single variant case (getFIX returns a vector if only one variant, screwing up further processing)
      if (!is.data.frame(fix_data)) {
        if (is.null(dim(fix_data))) {
          fix_data <- as.data.frame(t(fix_data), stringsAsFactors = FALSE)
        } else {
          fix_data <- as.data.frame(fix_data, stringsAsFactors = FALSE)
        }
      }
      
      gene_data$variants <- fix_data %>%
        mutate(
          POS = suppressWarnings(as.numeric(as.character(POS))),
          TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNP", "INDEL")
        ) %>% 
        filter(!is.na(POS))
      
      showNotification("Gene loaded successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
    
    
  })
  
  # Check if gene has been loaded
  
  output$gene_loaded <- reactive({
    return(!is.null(gene_data$variants))
  })
  outputOptions(output, "gene_loaded", suspendWhenHidden = FALSE)
  
  
  # Check if variants in gene
  
  output$has_variants <- reactive({
    req(gene_data$variants)
    return(nrow(gene_data$variants) > 0)
  })
  outputOptions(output, "has_variants", suspendWhenHidden = FALSE)
  
  # Gene location card
  
  output$gene_location_card <- renderUI({
    req(gene_data$location)
    
    
    chr <- gene_data$location$chr[1]
    start <- format(gene_data$location$start[1], big.mark = ",")
    end <- format(gene_data$location$end[1], big.mark = ",")
    length <- format(gene_data$location$end[1] - gene_data$location$start[1], big.mark = ",")
    
    div(
      style = "background-color: #f0f8ff; padding: 15px; border-radius: 5px; border-left: 5px solid #4CAF50;",
      h3(gene_data$current_gene, style = "margin-top: 0;"),
      p(strong("Location: "), sprintf("Chromosome %s: %s - %s bp", chr, start, end)),
      p(strong("Length: "), sprintf("%s bp", length), style = "margin-bottom: 0;")
    )
    
    
  })
  
  # Gene diagram with variants
  
  output$gene_diagram <- renderPlot({
    req(gene_data$features)
    
    
    variants_to_plot <- filtered_variants()
    
    plot_gene_structure(gene_data$features, variants_to_plot, gene_data$current_gene)
    
    
  })
  
  # Features table
  
  output$features_table <- renderDT({
    req(gene_data$features)
    
    
    gene_data$features %>%
      mutate(
        Length = end - start + 1,
        start = format(start, big.mark = ","),
        end = format(end, big.mark = ","),
        Length = format(Length, big.mark = ",")
      ) %>%
      datatable(
        options = list(
          pageLength = 10,
          dom = 't',
          ordering = FALSE
        ),
        rownames = FALSE,
        colnames = c("Feature", "Start", "End", "Length")
      )
    
    
  })
  
  # Homologs table
  
  output$homologs_table <- renderDT({
    req(gene_data$homologs)
    
    
    datatable(
      gene_data$homologs,
      options = list(
        pageLength = 5,
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE
    )
    
    
  })
  
  # Filtered variants based on user input
  
  filtered_variants <- reactive({
    req(gene_data$variants)
    
    
    variants <- gene_data$variants
    
    # Filter by variant type
    if (input$variant_type != "all") {
      variants <- variants %>% filter(TYPE == input$variant_type)
    }
    
    # Filter by coding regions if selected
    if (input$coding_only && !is.null(gene_data$features)) {
      coding_features <- gene_data$features %>%
        filter(feature_type %in% c("CDS", "exon"))
      
      variants <- variants %>%
        filter(sapply(POS, function(pos) {
          any(pos >= coding_features$start & pos <= coding_features$end)
        }))
    }
    
    return(variants)
    
    
  })
  
  # Filtered genotypes
  filtered_genotypes <- reactive({
    req(gene_data$vcf)
    
    variants <- filtered_variants()
    
    if (nrow(variants) == 0) return(NULL)
    
    # Extract genotypes
    gt_data <- extract.gt(gene_data$vcf, element = "GT")
    
    # Filter to matching variants
    variant_indices <- which(gene_data$variants$POS %in% variants$POS)
    gt_filtered <- gt_data[variant_indices, , drop = FALSE]
    
    
    # Combine with variant info
    result <- cbind(
      CHROM = variants$CHROM,
      POS = variants$POS,
      REF = variants$REF,
      ALT = variants$ALT,
      TYPE = variants$TYPE,
      gt_filtered
    )
    
    result <- as.data.frame(result)
    
    # Transpose the table so samples are rows
    result_t <- as.data.frame(t(result))
    colnames(result_t) <- paste0(result$CHROM, "_", result$POS)  # Variant positions as column names
    result_t <- result_t[-c(1:6), , drop = FALSE]  # Remove the metadata rows (CHROM, POS, REF, ALT, TYPE, MAF)
    result_t <- cbind(Sample = rownames(result_t), result_t)  # Add sample names as first column
    
    return(result_t)
  })
  
  # Variant summary
  output$variant_summary <- renderUI({
    req(gene_data$variants)
    
    all_vars <- gene_data$variants
    filtered_vars <- filtered_variants()
    
    total <- nrow(all_vars)
    filtered <- nrow(filtered_vars)
    snps <- sum(filtered_vars$TYPE == "SNP")
    indels <- sum(filtered_vars$TYPE == "INDEL")
    
    tagList(
      # Summary boxes
      fluidRow(
        column(4, 
               div(style = "text-align: center;",
                   h3(total, style = "color: #333; margin: 0;"),
                   p("Total Variants in Gene Region", style = "margin: 0;")
               )
        ),
        column(4,
               div(style = "text-align: center;",
                   h3(snps, style = "color: #2196F3; margin: 0;"),
                   p("SNPs", style = "margin: 0;")
               )
        ),
        column(4,
               div(style = "text-align: center;",
                   h3(indels, style = "color: #FF9800; margin: 0;"),
                   p("Indels", style = "margin: 0;")
               )
        )
      ),
      
      br(),
      
      # Variant details table
      h5("Variant Details"),
      DTOutput("variant_details_table")
    )
  })
  
  # Variant details table output
  output$variant_details_table <- renderDT({
    req(filtered_variants())
    
    variants <- filtered_variants() %>%
      select(CHROM, POS, REF, ALT, TYPE) %>%
      mutate(POS = format(POS, big.mark = ","))
    
    datatable(
      variants,
      options = list(
        pageLength = 10,
        dom = 't',
        scrollX = TRUE
      ),
      rownames = FALSE,
      colnames = c("Chromosome", "Position", "Reference", "Alternate", "Type")
    )
  })
  
  # Genotypes table
  output$genotypes_table <- renderDT({
    req(filtered_genotypes())
    
    datatable(
      filtered_genotypes(),
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      rownames = FALSE,
      filter = 'top'
    )
  })
  
  # Download genotypes as CSV
  
  output$download_genotypes <- downloadHandler(
    filename = function() {
      paste0(gene_data$current_gene, "*genotypes*", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_genotypes(), file, row.names = FALSE)
    }
  )
  
  # Download filtered VCF
  
  output$download_vcf <- downloadHandler(
    filename = function(file) {
      paste0(gene_data$current_gene, "*filtered*", Sys.Date(), ".vcf")
    },
    content = function(file) {
      req(gene_data$vcf)
      
      variants <- filtered_variants()
      variant_indices <- which(gene_data$variants$POS %in% variants$POS)
      
      vcf_filtered <- gene_data$vcf[variant_indices, ]
      
      # Writing vcf file as plain text because write.vcf() doesn't
      con <- file(file, open = "wt")
      
      writeLines(vcf_filtered@meta, con)
      writeLines(paste0("#", paste(colnames(vcf_filtered@fix), colnames(vcf_filtered@gt), sep = "\t", collapse = "\t")), con)
      write.table(cbind(vcf_filtered@fix, vcf_filtered@gt),
                  file = con,
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE)
      close(con)
    }
    
    
  )
}


# ===== Run the app =====
shinyApp(ui = ui, server = server)