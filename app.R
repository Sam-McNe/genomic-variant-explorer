# Variant Mining Tool by Sam McNeill | USDA-ARS | June 2026

# ========================================================================
# Project R library
# ========================================================================

PROJECT_R_LIB <- "/project/sugarbeet/samuel/Rlibs/R/RStudio-x86_64-pc-linux-gnu-library/4.4"

if (!dir.exists(PROJECT_R_LIB)) {
  stop(
    "Project R library not found: ",
    PROJECT_R_LIB,
    "\nThis app requires the shared project R library."
  )
}

.libPaths(unique(c(PROJECT_R_LIB, .libPaths())))

message("Using R library paths:")
message(paste(.libPaths(), collapse = "\n"))

# ========================================================================
# Packages
# ========================================================================

cran_packages <- c(
  "shiny",
  "ggplot2",
  "dplyr",
  "DT",
  "vcfR"
)

bioc_packages <- c(
  "Biostrings",
  "Rsamtools",
  "GenomicRanges",
  "rtracklayer",
  "msa"
)

all_packages <- c(cran_packages, bioc_packages)

missing_packages <- all_packages[
  !vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nExpected package library: ",
    PROJECT_R_LIB,
    "\nPlease install missing packages into the shared project library before running the app."
  )
}

invisible(lapply(all_packages, library, character.only = TRUE))


# ===== Directory Config =====

DATA_DIR        <- "/project/sugarbeet/samuel/VMT_Testing"
GFF_FILE        <- "./Bvulgarisssp_vulgaris_782_EL10.2_2.gene_exons.gff3"
VCF_FILE        <- "./mitchseq_maf1pct_maxmissing95pct_minMeanDP10_maxMeanDP500.vcf.gz"
ANNOTATION_FILE <- "./Bvulgarisssp_vulgaris_782_EL10.2_2.annotation_info.txt"
FASTA_FILE      <- "./Bvulgarisssp_vulgaris_782_EL10.2.fa"


# ===== Load GFF once at startup =====

message("Loading GFF file...")
GFF_DATA <- rtracklayer::import(GFF_FILE)
message("GFF loaded - ", length(GFF_DATA), " features")


# ===== Helper: Safe SnpEff ANN parser =====

ann_field <- function(ann, index) {
  if (is.na(ann) || !nzchar(ann)) return(NA_character_)
  
  first_ann <- strsplit(ann, ",", fixed = TRUE)[[1]][1]
  fields <- strsplit(first_ann, "|", fixed = TRUE)[[1]]
  
  if (length(fields) < index) return(NA_character_)
  fields[index]
}


# ===== Gene Diagram Function =====

plot_gene_structure <- function(features_df, variants_df = NULL, gene_name = "") {
  
  feature_hierarchy <- data.frame(
    feature_type = c(
      "gene", "mRNA", "CDS", "five_prime_UTR", "three_prime_UTR",
      "exon", "intron"
    ),
    track = c(1, 2, 3, 4, 5, 6, 7),
    color = c(
      "gray60", "gray40", "darkblue", "lightgreen", "lightgreen",
      "steelblue", "gray80"
    ),
    stringsAsFactors = FALSE
  )
  
  plot_data <- features_df %>%
    left_join(feature_hierarchy, by = "feature_type") %>%
    mutate(
      track = ifelse(is.na(track), max(feature_hierarchy$track) + 1, track),
      color = ifelse(is.na(color), "gray50", color)
    )
  
  gene_start <- min(plot_data$start)
  gene_end   <- max(plot_data$end)
  
  p <- ggplot(plot_data) +
    geom_rect(
      aes(
        xmin = start, xmax = end,
        ymin = track - 0.35, ymax = track + 0.35,
        fill = feature_type
      ),
      color = "black", linewidth = 0.3
    ) +
    geom_text(
      aes(
        x = gene_start - (gene_end - gene_start) * 0.02,
        y = track,
        label = feature_type
      ),
      hjust = 1, size = 3.5, fontface = "bold"
    ) +
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
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position    = "none",
      plot.title         = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x       = element_text(size = 11),
      axis.text.x        = element_text(size = 10)
    )
  
  if (!is.null(variants_df) && nrow(variants_df) > 0) {
    p <- p +
      geom_segment(
        data = variants_df,
        aes(
          x = POS, xend = POS,
          y = 0.5, yend = max(plot_data$track) + 0.5,
          color = TYPE
        ),
        alpha = 0.7, linewidth = 1
      ) +
      scale_color_manual(
        values = c("SNP" = "red", "INDEL" = "orange", "OTHER" = "gray"),
        name   = "Variant Type"
      )
    
    lof_variants <- variants_df[!is.na(variants_df$LOF) & variants_df$LOF != "", ]
    
    if (nrow(lof_variants) > 0) {
      p <- p +
        geom_point(
          data = lof_variants,
          aes(x = POS, y = max(plot_data$track) + 0.5),
          shape = 25, size = 5,
          fill = "#8B0000", color = "#8B0000",
          inherit.aes = FALSE
        ) +
        geom_text(
          data = lof_variants,
          aes(x = POS, y = max(plot_data$track) + 0.9),
          label = "LOF", size = 3, color = "#8B0000",
          fontface = "bold", inherit.aes = FALSE
        )
    }
  }
  
  return(p)
}


# ===== Build per-variant protein sequences =====

build_variant_proteins <- function(features_df, variants_df, gene_name,
                                   fasta_file, chr, gff_data) {
  
  fai_file <- paste0(fasta_file, ".fai")
  if (!file.exists(fai_file)) {
    Rsamtools::indexFa(fasta_file)
  }
  
  fa <- Rsamtools::FaFile(fasta_file)
  open(fa)
  on.exit(close(fa), add = TRUE)
  
  gene_row <- gff_data[
    gff_data$type == "gene" &
      grepl(gene_name, as.character(gff_data$ID), fixed = TRUE)
  ]
  
  if (length(gene_row) == 0) {
    warning("Gene not found in GFF")
    return(NULL)
  }
  
  strand <- as.character(GenomicRanges::strand(gene_row))[1]
  
  cds_features <- features_df %>%
    filter(feature_type == "CDS") %>%
    distinct(start, end, .keep_all = TRUE) %>%
    {
      if (strand == "+") arrange(., start) else arrange(., desc(start))
    }
  
  if (nrow(cds_features) == 0) {
    warning("No CDS features for gene")
    return(NULL)
  }
  
  ref_exons <- lapply(seq_len(nrow(cds_features)), function(i) {
    gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(
        start = cds_features$start[i],
        end = cds_features$end[i]
      )
    )
    
    raw <- Rsamtools::scanFa(fa, gr)[[1]]
    
    if (strand == "-") {
      Biostrings::reverseComplement(raw)
    } else {
      raw
    }
  })
  
  ref_cds_dna <- do.call(Biostrings::xscat, ref_exons)
  
  n <- nchar(ref_cds_dna)
  if (n < 3) {
    warning("Reference CDS is shorter than one codon")
    return(NULL)
  }
  
  ref_cds_trimmed <- Biostrings::subseq(ref_cds_dna, 1, n - n %% 3)
  ref_protein <- Biostrings::translate(ref_cds_trimmed, if.fuzzy.codon = "solve")
  
  protein_list <- list()
  protein_list[["Reference"]] <- ref_protein
  
  message(">>> Reference protein length: ", nchar(ref_protein))
  
  cds_ranges <- features_df %>%
    filter(feature_type == "CDS")
  
  coding_variants <- variants_df %>%
    filter(
      sapply(POS, function(pos) {
        any(pos >= cds_ranges$start & pos <= cds_ranges$end)
      }),
      !grepl(",", ALT, fixed = TRUE),
      grepl("^[ACGTNacgtn]+$", REF),
      grepl("^[ACGTNacgtn]+$", ALT)
    )
  
  if (nrow(coding_variants) == 0) {
    return(Biostrings::AAStringSet(protein_list))
  }
  
  message("Found ", nrow(coding_variants), " coding variants. Building alternate proteins...")
  
  cds_features$length <- cds_features$end - cds_features$start + 1
  cds_features$cumulative_len <- cumsum(c(0, cds_features$length[-nrow(cds_features)]))
  
  ref_cds_string <- as.character(ref_cds_dna)
  
  for (i in seq_len(nrow(coding_variants))) {
    var <- coding_variants[i, ]
    
    var_pos <- as.integer(var$POS)
    
    exon_idx <- which(var_pos >= cds_features$start & var_pos <= cds_features$end)
    if (length(exon_idx) == 0) next
    
    current_exon <- cds_features[exon_idx[1], ]
    
    pos_in_exon <- if (strand == "+") {
      var_pos - current_exon$start + 1
    } else {
      current_exon$end - var_pos + 1
    }
    
    pos_in_cds <- current_exon$cumulative_len + pos_in_exon
    
    ref_allele_sense <- if (strand == "-") {
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(var$REF)))
    } else {
      toupper(var$REF)
    }
    
    alt_allele_sense <- if (strand == "-") {
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(var$ALT)))
    } else {
      toupper(var$ALT)
    }
    
    fasta_ref <- substr(
      ref_cds_string,
      pos_in_cds,
      pos_in_cds + nchar(ref_allele_sense) - 1
    )
    
    if (toupper(fasta_ref) != toupper(ref_allele_sense)) {
      warning(
        "REF mismatch at POS ", var_pos,
        ". FASTA has ", fasta_ref,
        ", VCF has ", ref_allele_sense,
        ". Skipping."
      )
      next
    }
    
    alt_cds_string <- paste0(
      substr(ref_cds_string, 1, pos_in_cds - 1),
      alt_allele_sense,
      substr(
        ref_cds_string,
        pos_in_cds + nchar(ref_allele_sense),
        nchar(ref_cds_string)
      )
    )
    
    alt_cds_dna <- Biostrings::DNAString(alt_cds_string)
    
    n_alt <- nchar(alt_cds_dna)
    if (n_alt < 3) next
    
    alt_cds_trimmed <- Biostrings::subseq(alt_cds_dna, 1, n_alt - n_alt %% 3)
    alt_protein <- Biostrings::translate(alt_cds_trimmed, if.fuzzy.codon = "solve")
    
    hgvs_p <- var[["HGVS.p"]]
    
    var_name <- if (!is.na(hgvs_p) && nzchar(hgvs_p)) {
      hgvs_p
    } else {
      paste0("var_", var$POS, "_", var$REF, "_", var$ALT)
    }
    
    var_name <- gsub("[^A-Za-z0-9_.-]", "_", var_name)
    
    message(
      ">>> Variant ", var_pos,
      " (", var_name, ") protein length: ",
      nchar(alt_protein)
    )
    
    protein_list[[var_name]] <- alt_protein
  }
  
  seq_chars <- vapply(protein_list, as.character, character(1))
  keep <- !duplicated(seq_chars)
  
  seq_chars <- seq_chars[keep]
  names(seq_chars) <- make.unique(names(seq_chars))
  
  protein_set <- Biostrings::AAStringSet(seq_chars)
  
  message("Returning ", length(protein_set), " unique protein sequences for alignment.")
  
  protein_set
}


# ========================================================================
# Alignment viewer helpers
# ========================================================================

split_aa <- function(x) {
  strsplit(as.character(x), split = "", fixed = TRUE)[[1]]
}

mask_after_first_stop <- function(chars, stop_symbol = "*", mask_symbol = "-") {
  stop_pos <- which(chars == stop_symbol)
  
  if (length(stop_pos) == 0) {
    return(chars)
  }
  
  first_stop <- stop_pos[1]
  
  if (first_stop < length(chars)) {
    chars[(first_stop + 1):length(chars)] <- mask_symbol
  }
  
  chars
}

aa_class_for_display <- function(
    display_char,
    ref_char,
    raw_char,
    after_stop = FALSE,
    is_ref = FALSE
) {
  if (is_ref) {
    return("aa-cell aa-ref")
  }
  
  if (after_stop) {
    return("aa-cell aa-after-stop")
  }
  
  if (raw_char == "*") {
    return("aa-cell aa-stop")
  }
  
  if (raw_char == "-") {
    return("aa-cell aa-gap")
  }
  
  if (display_char == "." && raw_char == ref_char) {
    return("aa-cell aa-identical")
  }
  
  if (raw_char != ref_char) {
    return("aa-cell aa-substitution")
  }
  
  "aa-cell"
}


# ========================================================================
# Align proteins with Clustal Omega and build tree
# ========================================================================

build_msa <- function(protein_set) {
  if (is.null(protein_set) || length(protein_set) < 2) {
    warning("Need at least 2 unique protein sequences to build an alignment.")
    return(NULL)
  }
  
  message("Aligning protein sequences with Clustal Omega...")
  
  aln <- msa::msa(
    protein_set,
    method = "ClustalOmega",
    order = "input"
  )
  
  aln_set <- as(aln, "AAStringSet")
  
  # Preserve names from original protein set where possible
  if (length(names(protein_set)) == length(aln_set)) {
    names(aln_set) <- names(protein_set)
  }
  
  # Ensure every sequence has a safe unique name
  bad_names <- is.na(names(aln_set)) | names(aln_set) == ""
  
  if (any(bad_names)) {
    names(aln_set)[bad_names] <- paste0("seq_", which(bad_names))
  }
  
  names(aln_set) <- make.unique(
    gsub("[^A-Za-z0-9_.-]", "_", names(aln_set))
  )
  
  list(
    alignment_set = aln_set
  )
}


# ===== UI =====

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .bioedit-scroll-container {
        overflow-x: auto;
        overflow-y: auto;
        max-width: 100%;
        max-height: 700px;
        border: 1px solid #ddd;
        background-color: white;
        padding: 8px;
        font-family: 'Courier New', monospace;
        font-size: 13px;
        white-space: nowrap;
      }

      .aa-ref {
        color: #111;
        background-color: #f7f7f7;
      }

      .aa-identical {
        color: #888;
        background-color: white;
      }

      .aa-substitution {
        color: #000;
        background-color: #ffd966;
        font-weight: bold;
      }

      .aa-gap {
        color: #999;
        background-color: #eeeeee;
      }

      .aa-stop {
        color: white;
        background-color: #cc0000;
        font-weight: bold;
      }
      
      .aa-after-stop {
        color: #777;
        background-color: #f1f1f1;
      }

      .bioedit-position-ruler {
        position: sticky;
        top: 0;
        z-index: 5;
        background-color: white;
        white-space: nowrap;
        color: #666;
        font-family: 'Courier New', monospace;
        font-size: 11px;
        margin-bottom: 4px;
        border-bottom: 1px solid #ddd;
      }
      
      .position-ruler-seq {
        font-size: 0;
        letter-spacing: 0;
        word-spacing: 0;
      }

      .position-block {
        display: inline-block;
        height: 16px;
        line-height: 16px;
        text-align: right;
        color: #337ab7;
        font-weight: bold;
        font-size: 11px;
        box-sizing: border-box;
        padding-right: 1px;
        vertical-align: top;
        letter-spacing: normal;
        word-spacing: normal;
      }

      .bioedit-row {
        display: flex;
        align-items: stretch;
        white-space: nowrap;
        line-height: 1.6em;
        min-width: max-content;
        width: max-content;
        position: relative;
      }
      
      .bioedit-label {
        display: block;
        position: sticky;
        left: 0;
        z-index: 200;
      
        flex: 0 0 220px;
        width: 220px;
        min-width: 220px;
        max-width: 220px;
      
        overflow: hidden;
        text-overflow: ellipsis;
        vertical-align: top;
        font-weight: bold;
        color: #333;
        white-space: nowrap;
      
        background-color: #ffffff;
        border-right: 1px solid #cccccc;
        padding-left: 6px;
        padding-right: 6px;
        margin-right: 0;
        box-sizing: border-box;
        box-shadow: 12px 0 0 0 #ffffff;
      }
      
      .bioedit-ref-label {
        color: #000000;
        background-color: #f0f0f0;
        z-index: 220;
        box-shadow: 12px 0 0 0 #f0f0f0;
      }
      
      .bioedit-position-ruler .bioedit-label {
        background-color: #ffffff;
        color: #337ab7;
        z-index: 230;
        box-shadow: 12px 0 0 0 #ffffff;
      }
      
      .bioedit-seq {
        display: block;
        flex: 0 0 auto;
        white-space: nowrap;
        vertical-align: top;
        min-width: max-content;
      
        font-size: 0;
        letter-spacing: 0;
        word-spacing: 0;
      
        position: relative;
        z-index: 1;
        margin-left: 0;
      }
      
      .aa-cell {
        display: inline-block;
        width: 14px;
        min-width: 14px;
        max-width: 14px;
        text-align: center;
        margin: 0;
        padding: 0;
        border-radius: 2px;
        white-space: pre;
        box-sizing: border-box;
      
        font-size: 13px;
        letter-spacing: normal;
        word-spacing: normal;
        vertical-align: top;
      
        position: relative;
        z-index: 1;
      }
            
      
      
      
    "))
  ),
  
  titlePanel("Variant Mining Tool"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      textInput(
        "gene_name",
        "Enter Gene Name:",
        placeholder = "e.g., Bevul.2G077600"
      ),
      
      actionButton(
        "submit",
        "Load Gene",
        class = "btn-primary",
        style = "width: 100%; margin-bottom: 20px;"
      ),
      
      hr(),
      
      h5("Variant Filters"),
      
      selectInput(
        "variant_type",
        "Variant Type:",
        choices = c("All" = "all", "SNP" = "SNP", "INDEL" = "INDEL"),
        selected = "all"
      ),
      
      checkboxInput("coding_only", "Coding variants only", value = TRUE),
      
      hr(),
      
      downloadButton(
        "download_genotypes",
        "Download Genotypes",
        style = "width: 100%; margin-bottom: 10px;"
      ),
      
      downloadButton(
        "download_vcf",
        "Download VCF",
        style = "width: 100%; margin-bottom: 10px;"
      ),
      
      hr(),
      
      h5("Top Homologs"),
      uiOutput("homologs_sidebar")
    ),
    
    mainPanel(
      width = 9,
      
      uiOutput("gene_location_card"),
      br(),
      
      tabsetPanel(
        tabPanel(
          "Gene & Variants",
          br(),
          plotOutput("gene_diagram", height = "350px"),
          br(),
          wellPanel(
            h4("Variant Summary"),
            uiOutput("variant_summary")
          )
        ),
        
        tabPanel(
          "Genotypes",
          br(),
          
          textInput(
            "genotype_col_search",
            "Search variant name / position:",
            placeholder = "e.g. 17160256"
          ),
          
          p(
            "Note: do not include commas in the position number.",
            style = "color: gray; font-size: 12px;"
          ),
          
          conditionalPanel(
            condition = "output.gene_loaded && output.has_variants",
            DTOutput("genotypes_table")
          ),
          
          conditionalPanel(
            condition = "output.gene_loaded && !output.has_variants",
            div(
              style = "text-align: center; padding: 40px; color: #999;",
              h4("No variants found in this gene", style = "color: #FF9800;"),
              p("Select a different gene.")
            )
          )
        ),
        
        tabPanel(
          "Protein Alignment",
          
          h4("Protein Variant Difference Viewer"),
          
          fluidRow(
            column(
              3,
              checkboxInput(
                "bioedit_diff_mode",
                "Show only variant AA changes",
                value = TRUE
              )
            )
          ),
          
          uiOutput("bioedit_alignment_stats"),
          
          downloadButton(
            "download_aln_fasta",
            "Download aligned protein FASTA",
            style = "margin-bottom: 10px;"
          ),
          
          div(
            class = "bioedit-scroll-container",
            uiOutput("bioedit_alignment_viewer")
          )
        )
    )
  )
)
)

# ===== Server =====

server <- function(input, output, session) {
  
  gene_data <- reactiveValues(
    location = NULL,
    features = NULL,
    homologs = NULL,
    vcf = NULL,
    variants = NULL,
    current_gene = NULL,
    msa_result = NULL,
  )
  
  
  observeEvent(input$submit, {
    req(input$gene_name)
    
    # Reset all reactive values
    lapply(names(gene_data), function(x) gene_data[[x]] <- NULL)
    
    gene_data$current_gene <- input$gene_name
    
    tryCatch({
      
      showNotification(
        "Gathering data...",
        type = "message",
        duration = NULL,
        id = "generating"
      )
      
      bash_script_path <- "/project/sugarbeet/samuel/VMT_Testing/variant_mining_tool_2.sh"
      
      result <- system2(
        "bash",
        args = c(bash_script_path, gene_data$current_gene),
        stdout = TRUE,
        stderr = TRUE,
        wait = TRUE
      )
      
      exit_code <- attr(result, "status")
      
      if (!is.null(exit_code) && exit_code != 0) {
        removeNotification(id = "generating")
        showNotification(
          "Gene not found. Please check the name and try again.",
          type = "error",
          duration = 5
        )
        return()
      }
      
      removeNotification(id = "generating")
      
      location_file <- file.path(
        DATA_DIR,
        "output",
        paste0(gene_data$current_gene, ".region.txt")
      )
      
      if (!file.exists(location_file)) {
        stop("Data files could not be found.")
      }
      
      gene_data$location <- read.table(
        location_file,
        header = FALSE,
        col.names = c("chr", "start", "end")
      )
      
      features_file <- file.path(
        DATA_DIR,
        "output",
        paste0(gene_data$current_gene, ".features.txt")
      )
      
      gene_data$features <- read.table(
        features_file,
        header = FALSE,
        col.names = c("feature_type", "start", "end")
      )
      
      homologs_file <- file.path(
        DATA_DIR,
        "output",
        paste0(gene_data$current_gene, ".annotation.txt")
      )
      
      gene_data$homologs <- read.table(
        homologs_file,
        header = TRUE,
        sep = "\t",
        quote = ""
      )
      
      colnames(gene_data$homologs) <- c(
        "Best hit Arabidopsis name",
        "Best hit Arabidopsis definition",
        "Best hit Clamydomonas name",
        "Best hit Clamydomonas definition",
        "Best hit Rice name",
        "Best hit Rice definition"
      )
      
      vcf_file <- file.path(
        DATA_DIR,
        "output",
        paste0(gene_data$current_gene, ".vcf")
      )
      
      variant_lines <- length(grep(
        "^#",
        readLines(vcf_file),
        invert = TRUE,
        value = TRUE
      ))
      
      if (variant_lines == 0) {
        
        gene_data$vcf <- NULL
        gene_data$variants <- data.frame(
          CHROM = character(),
          POS = numeric(),
          ID = character(),
          REF = character(),
          ALT = character(),
          QUAL = character(),
          FILTER = character(),
          INFO = character(),
          TYPE = character(),
          AF = numeric(),
          LOF = character(),
          Annotation = character(),
          Annotation_Impact = character(),
          `HGVS.c` = character(),
          `HGVS.p` = character(),
          check.names = FALSE
        )
        
        showNotification(
          "Gene loaded. No variants found.",
          type = "warning",
          duration = 4
        )
        
      } else {
        
        gene_data$vcf <- read.vcfR(vcf_file, verbose = FALSE)
        
        fix_data <- as.data.frame(
          getFIX(gene_data$vcf, getINFO = TRUE),
          stringsAsFactors = FALSE
        )
        
        gene_data$variants <- fix_data %>%
          mutate(
            POS = as.numeric(POS),
            TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNP", "INDEL"),
            AF = suppressWarnings(as.numeric(
              vcfR::extract.info(gene_data$vcf, element = "AF")
            )),
            LOF = vcfR::extract.info(gene_data$vcf, element = "LOF"),
            ANN = vcfR::extract.info(gene_data$vcf, element = "ANN"),
            Annotation = vapply(ANN, ann_field, character(1), index = 2),
            Annotation_Impact = vapply(ANN, ann_field, character(1), index = 3),
            `HGVS.c` = vapply(ANN, ann_field, character(1), index = 10),
            `HGVS.p` = vapply(ANN, ann_field, character(1), index = 11)
          ) %>%
          filter(!is.na(POS))
        
        showNotification(
          "Building variant-specific proteins...",
          type = "message",
          duration = NULL,
          id = "msa_msg"
        )
        
        protein_set <- build_variant_proteins(
          features_df = gene_data$features,
          variants_df = gene_data$variants,
          gene_name = gene_data$current_gene,
          fasta_file = FASTA_FILE,
          chr = gene_data$location$chr[1],
          gff_data = GFF_DATA
        )
        
        removeNotification(id = "msa_msg")
        
        if (!is.null(protein_set) && length(protein_set) > 1) {
          
          showNotification(
            "Aligning proteins and building tree...",
            type = "message",
            duration = NULL,
            id = "aln_msg"
          )
          
          gene_data$msa_result <- build_msa(protein_set)
          
          removeNotification(id = "aln_msg")
          
          if (!is.null(gene_data$msa_result)) {
            showNotification("Alignment complete!", type = "message", duration = 4)
          } else {
            showNotification(
              "Protein sequences were generated, but alignment failed.",
              type = "warning",
              duration = 6
            )
          }
          
        } else {
          showNotification(
            "No coding variants found to generate an alignment.",
            type = "warning",
            duration = 5
          )
        }
        
        showNotification(
          "Gene loaded successfully!",
          type = "message",
          duration = 3
        )
      }
      
    }, error = function(e) {
      
      lapply(c("generating", "msa_msg", "aln_msg"), removeNotification)
      
      showNotification(
        paste("Error:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  
  output$gene_loaded <- reactive({
    !is.null(gene_data$variants)
  })
  outputOptions(output, "gene_loaded", suspendWhenHidden = FALSE)
  
  
  output$has_variants <- reactive({
    req(gene_data$variants)
    nrow(gene_data$variants) > 0
  })
  outputOptions(output, "has_variants", suspendWhenHidden = FALSE)
  
  output$gene_location_card <- renderUI({
    req(gene_data$location)
    
    div(
      style = "background-color: #f0f8ff; padding: 15px; border-radius: 5px; border-left: 5px solid #4CAF50;",
      
      h3(gene_data$current_gene, style = "margin-top: 0;"),
      
      p(
        strong("Location: "),
        sprintf(
          "Chromosome %s: %s - %s bp",
          gene_data$location$chr[1],
          format(gene_data$location$start[1], big.mark = ","),
          format(gene_data$location$end[1], big.mark = ",")
        )
      ),
      
      p(
        strong("Length: "),
        sprintf(
          "%s bp",
          format(
            gene_data$location$end[1] - gene_data$location$start[1],
            big.mark = ","
          )
        ),
        style = "margin-bottom: 0;"
      )
    )
  })
  
  
  output$gene_diagram <- renderPlot({
    req(gene_data$features)
    plot_gene_structure(
      gene_data$features,
      filtered_variants(),
      gene_data$current_gene
    )
  })
  
  
  output$homologs_sidebar <- renderUI({
    req(gene_data$homologs)
    
    h <- gene_data$homologs
    
    div(
      style = "font-size: 12px;",
      
      div(
        style = "margin-bottom: 8px;",
        strong("Arabidopsis: "),
        div(h[1, 1], style = "color: #555;"),
        div(h[1, 2], style = "color: #777; font-style: italic;")
      ),
      
      div(
        style = "margin-bottom: 8px;",
        strong("Chlamydomonas: "),
        div(h[1, 3], style = "color: #555;"),
        div(h[1, 4], style = "color: #777; font-style: italic;")
      ),
      
      div(
        style = "margin-bottom: 8px;",
        strong("Rice: "),
        div(h[1, 5], style = "color: #555;"),
        div(h[1, 6], style = "color: #777; font-style: italic;")
      )
    )
  })
  
  
  filtered_variants <- reactive({
    req(gene_data$variants)
    
    variants <- gene_data$variants
    
    if (nrow(variants) == 0) {
      return(variants)
    }
    
    if (input$variant_type != "all") {
      variants <- variants %>%
        filter(TYPE == input$variant_type)
    }
    
    if (input$coding_only) {
      req(gene_data$features)
      
      coding_features <- gene_data$features %>%
        filter(feature_type == "CDS")
      
      if (nrow(coding_features) > 0) {
        variants <- variants %>%
          filter(
            sapply(POS, function(pos) {
              any(pos >= coding_features$start & pos <= coding_features$end)
            })
          )
      }
    }
    
    variants
  })
  
  
  output$variant_summary <- renderUI({
    req(filtered_variants())
    
    all_vars <- gene_data$variants
    filtered_vars <- filtered_variants()
    
    total <- nrow(filtered_vars)
    snps <- sum(filtered_vars$TYPE == "SNP", na.rm = TRUE)
    indels <- sum(filtered_vars$TYPE == "INDEL", na.rm = TRUE)
    
    tagList(
      fluidRow(
        column(
          4,
          div(
            style = "text-align: center;",
            h3(total, style = "color: #333; margin: 0;"),
            p("Total Variants", style = "margin: 0;")
          )
        ),
        
        column(
          4,
          div(
            style = "text-align: center;",
            h3(snps, style = "color: #2196F3; margin: 0;"),
            p("SNPs", style = "margin: 0;")
          )
        ),
        
        column(
          4,
          div(
            style = "text-align: center;",
            h3(indels, style = "color: #FF9800; margin: 0;"),
            p("Indels", style = "margin: 0;")
          )
        )
      ),
      
      br(),
      
      DTOutput("variant_details_table")
    )
  })
  
  
  output$variant_details_table <- renderDT({
    req(filtered_variants())
    
    variants <- filtered_variants()
    
    if (nrow(variants) == 0) {
      return(datatable(
        data.frame(Message = "No variants found"),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    variants_for_table <- variants %>%
      mutate(
        POS_display = format(POS, big.mark = ","),
        LOF_flag = ifelse(is.na(LOF) | LOF == "", "", "⚠ LOF"),
        Impact_order = case_when(
          Annotation_Impact == "HIGH" ~ 1,
          Annotation_Impact == "MODERATE" ~ 2,
          Annotation_Impact == "LOW" ~ 3,
          Annotation_Impact == "MODIFIER" ~ 4,
          TRUE ~ 5
        )
      ) %>%
      select(
        POS_display, REF, ALT, TYPE, AF,
        Annotation_Impact, Annotation, LOF_flag, Impact_order
      )
    
    datatable(
      variants_for_table,
      caption = "Double-click any row for full variant details.",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        scrollY = "400px",
        scrollCollapse = TRUE,
        dom = "frtip",
        columnDefs = list(
          list(visible = FALSE, targets = 8),
          list(orderData = 8, targets = 5)
        ),
        initComplete = JS(
          "function(settings, json) {",
          "  var table = this;",
          "  $(this.api().table().body()).on('dblclick', 'tr', function() {",
          "    var rowData = table.api().row(this).data();",
          "    Shiny.setInputValue('variant_dblclick', rowData[0], {priority: 'event'});",
          "  });",
          "}"
        )
      ),
      rownames = FALSE,
      colnames = c(
        "Position", "REF", "ALT", "Type", "Allele Freq",
        "Impact", "Annotation", "Loss Of Function", ""
      )
    )
  })
  
  
  observeEvent(input$variant_dblclick, {
    req(input$variant_dblclick)
    
    pos_clicked <- gsub(",", "", input$variant_dblclick)
    
    variant <- filtered_variants() %>%
      filter(as.character(POS) == pos_clicked) %>%
      dplyr::slice(1)
    
    req(nrow(variant) > 0)
    
    showModal(modalDialog(
      title = paste0(
        "Variant Detail — ",
        variant$CHROM,
        " : ",
        format(as.numeric(variant$POS), big.mark = ",")
      ),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      
      fluidRow(
        column(
          6,
          h5(
            "Basic Info",
            style = "font-weight:bold; border-bottom:1px solid #ddd; padding-bottom:5px;"
          ),
          tags$table(
            class = "table table-sm table-bordered",
            tags$tr(tags$td(strong("Chromosome")), tags$td(variant$CHROM)),
            tags$tr(tags$td(strong("Position")), tags$td(format(as.numeric(variant$POS), big.mark = ","))),
            tags$tr(tags$td(strong("REF")), tags$td(variant$REF)),
            tags$tr(tags$td(strong("ALT")), tags$td(variant$ALT)),
            tags$tr(tags$td(strong("Type")), tags$td(variant$TYPE)),
            tags$tr(tags$td(strong("Allele Freq")), tags$td(round(variant$AF, 3)))
          )
        ),
        
        column(
          6,
          h5(
            "Functional Annotation",
            style = "font-weight:bold; border-bottom:1px solid #ddd; padding-bottom:5px;"
          ),
          tags$table(
            class = "table table-sm table-bordered",
            tags$tr(tags$td(strong("Annotation")), tags$td(variant$Annotation)),
            tags$tr(
              tags$td(strong("Impact")),
              tags$td(
                span(
                  variant$Annotation_Impact,
                  style = switch(
                    variant$Annotation_Impact,
                    "HIGH"     = "color:red; font-weight:bold;",
                    "MODERATE" = "color:orange; font-weight:bold;",
                    "LOW"      = "color:green;",
                    "MODIFIER" = "color:gray;",
                    ""
                  )
                )
              )
            ),
            tags$tr(tags$td(strong("HGVS.c")), tags$td(variant[["HGVS.c"]])),
            tags$tr(tags$td(strong("HGVS.p")), tags$td(variant[["HGVS.p"]])),
            tags$tr(
              tags$td(strong("LOF")),
              tags$td(ifelse(is.na(variant$LOF) || variant$LOF == "", "None", variant$LOF))
            )
          )
        )
      ),
      
      hr(),
      
      h5("All Annotations", style = "font-weight: bold;"),
      p(
        "All SnpEff annotations for this variant.",
        style = "color: gray; font-size: 12px;"
      ),
      
      DTOutput("modal_ann_table")
    ))
  })
  
  
  output$modal_ann_table <- renderDT({
    req(input$variant_dblclick)
    
    pos_clicked <- gsub(",", "", input$variant_dblclick)
    
    variant <- filtered_variants() %>%
      filter(as.character(POS) == pos_clicked) %>%
      dplyr::slice(1)
    
    req(nrow(variant) > 0)
    
    ann_field_value <- gene_data$variants %>%
      filter(POS == variant$POS) %>%
      pull(INFO) %>%
      dplyr::first()
    
    if (is.na(ann_field_value) || !grepl("ANN=", ann_field_value)) {
      return(datatable(
        data.frame(Message = "No ANN field found for this variant."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    ann_field_value <- sub(".*ANN=([^;]+).*", "\\1", ann_field_value)
    ann_entries <- strsplit(ann_field_value, ",", fixed = TRUE)[[1]]
    
    ann_df <- do.call(rbind, lapply(ann_entries, function(entry) {
      fields <- strsplit(entry, "|", fixed = TRUE)[[1]]
      length(fields) <- 16
      as.data.frame(t(fields), stringsAsFactors = FALSE)
    })) %>%
      setNames(c(
        "Allele", "Annotation", "Impact", "Gene_Name", "Gene_ID",
        "Feature_Type", "Feature_ID", "BioType", "Rank", "HGVS.c",
        "HGVS.p", "cDNA", "CDS", "AA", "Distance", "Errors"
      )) %>%
      select(
        Annotation, Impact, Gene_Name, Feature_ID,
        BioType, Rank, `HGVS.c`, `HGVS.p`
      )
    
    datatable(
      ann_df,
      options = list(dom = "t", pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  
  filtered_genotypes <- reactive({
    req(gene_data$vcf, filtered_variants())
    
    variants <- filtered_variants()
    
    if (nrow(variants) == 0) {
      return(NULL)
    }
    
    gt_data <- vcfR::extract.gt(gene_data$vcf, element = "GT")
    vcf_pos <- vcfR::getFIX(gene_data$vcf)[, "POS"]
    
    variant_indices <- which(as.numeric(vcf_pos) %in% variants$POS)
    
    if (length(variant_indices) == 0) {
      return(NULL)
    }
    
    gt_filtered <- gt_data[variant_indices, , drop = FALSE]
    
    info_filtered <- as.data.frame(
      vcfR::getFIX(gene_data$vcf)[variant_indices, c("CHROM", "POS", "REF", "ALT")],
      stringsAsFactors = FALSE
    )
    
    info_filtered$TYPE <- variants$TYPE[
      match(as.numeric(info_filtered$POS), variants$POS)
    ]
    
    result <- cbind(info_filtered, gt_filtered)
    
    result_t <- as.data.frame(t(result))
    
    colnames(result_t) <- paste0(result$CHROM, "_", result$POS)
    
    result_t <- result_t[-c(1:5), , drop = FALSE]
    
    result_t <- cbind(
      Sample = rownames(result_t),
      result_t
    )
    
    rownames(result_t) <- NULL
    
    result_t
  })
  
  output$genotypes_table <- renderDT({
    gt <- filtered_genotypes()
    
    if (is.null(gt)) {
      return(datatable(
        data.frame(Message = "No genotypes available for the selected variants."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    if (!is.null(input$genotype_col_search) && nzchar(input$genotype_col_search)) {
      matching_cols <- c(
        "Sample",
        grep(input$genotype_col_search, colnames(gt), value = TRUE)
      )
      
      gt <- gt[, colnames(gt) %in% matching_cols, drop = FALSE]
    }
    
    datatable(
      gt,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = "rtip"
      ),
      rownames = FALSE,
      filter = "top"
    )
  })
  
  
  
  output$bioedit_alignment_stats <- renderUI({
    req(gene_data$msa_result)
    req(gene_data$msa_result$alignment_set)
    
    aln_set <- gene_data$msa_result$alignment_set
    
    seq_names <- names(aln_set)
    seqs <- as.character(aln_set)
    
    if (is.null(seq_names) || any(is.na(seq_names)) || any(seq_names == "")) {
      seq_names <- paste0("seq_", seq_along(seqs))
    }
    
    ref_idx <- which(seq_names == "Reference")
    
    if (length(ref_idx) == 0) {
      ref_idx <- 1
    } else {
      ref_idx <- ref_idx[1]
    }
    
    ref_seq <- seqs[ref_idx]
    
    differs_from_ref <- seqs != ref_seq
    n_diff <- sum(differs_from_ref[-ref_idx], na.rm = TRUE)
    
    n_total <- length(seqs)
    aln_len <- max(nchar(seqs))
    
    div(
      style = "background-color: #f8f9fa; border-left: 4px solid #337ab7; padding: 10px; margin-bottom: 10px;",
      strong("Protein alignment overview: "),
      
      paste0(
        n_total, " unique aligned protein sequences, including Reference; ",
        n_diff, " unique variant protein sequences differ from Reference; ",
        aln_len, " amino-acid alignment columns."
      ),
      
      br(),
      
      em("*Synonymous variants may be omitted because they produce the same protein sequence. "),
      em("Identical residues are shown as dots.")
    )
  })
  
  
  output$bioedit_alignment_viewer <- renderUI({
    req(gene_data$msa_result)
    req(gene_data$msa_result$alignment_set)
    
    aln_set <- gene_data$msa_result$alignment_set
    
    validate(
      need(length(aln_set) > 0, "No aligned protein sequences available.")
    )
    
    seq_names <- names(aln_set)
    seqs <- as.character(aln_set)
    
    if (is.null(seq_names) || any(is.na(seq_names)) || any(seq_names == "")) {
      seq_names <- paste0("seq_", seq_along(seqs))
    }
    
    ref_idx <- which(seq_names == "Reference")
    
    if (length(ref_idx) == 0) {
      ref_idx <- 1
    } else {
      ref_idx <- ref_idx[1]
    }
    
    ref_seq <- seqs[ref_idx]
    ref_chars <- split_aa(ref_seq)
    
    aln_len <- length(ref_chars)
    
    seq_chars_list <- lapply(seqs, function(s) {
      x <- split_aa(s)
      
      if (length(x) < aln_len) {
        x <- c(x, rep("-", aln_len - length(x)))
      }
      
      if (length(x) > aln_len) {
        x <- x[seq_len(aln_len)]
      }
      
      x
    })
    
    rows <- list()
    aa_cell_px <- 14
    
    ruler_blocks <- list()
    
    full_blocks <- floor(aln_len / 10)
    remainder <- aln_len %% 10
    
    if (full_blocks > 0) {
      for (b in seq_len(full_blocks)) {
        pos_label <- b * 10
        
        ruler_blocks[[length(ruler_blocks) + 1]] <- span(
          class = "position-block",
          style = paste0("width: ", 10 * aa_cell_px, "px;"),
          as.character(pos_label)
        )
      }
    }
    
    if (remainder > 0) {
      ruler_blocks[[length(ruler_blocks) + 1]] <- span(
        class = "position-block",
        style = paste0("width: ", remainder * aa_cell_px, "px;"),
        ""
      )
    }
    
    rows[[length(rows) + 1]] <- div(
      class = "bioedit-row bioedit-position-ruler",
      span(class = "bioedit-label", "Position"),
      span(
        class = "bioedit-seq position-ruler-seq",
        ruler_blocks
      )
    )
    
    for (i in seq_along(seqs)) {
      raw_chars <- seq_chars_list[[i]]
      is_ref <- i == ref_idx
      
      # Identify positions after the first stop codon.
      # This is always enabled for variant rows.
      stop_pos <- which(raw_chars == "*")
      after_stop <- rep(FALSE, aln_len)
      
      if (!is_ref && length(stop_pos) > 0 && stop_pos[1] < aln_len) {
        after_stop[(stop_pos[1] + 1):aln_len] <- TRUE
      }
      
      # Always mask downstream residues after first stop in variant rows.
      masked_chars <- raw_chars
      
      if (!is_ref) {
        masked_chars <- mask_after_first_stop(masked_chars, mask_symbol = "-")
      }
      
      display_chars <- masked_chars

      if (isTRUE(input$bioedit_diff_mode) && !is_ref) {
        same_as_ref <- raw_chars == ref_chars
        
        display_chars[same_as_ref] <- "."
        
        if (any(after_stop)) {
          display_chars[after_stop] <- "-"
        }
      }
      
      aa_spans <- lapply(seq_len(aln_len), function(pos) {
        raw_char <- raw_chars[pos]
        display_char <- display_chars[pos]
        
        css_class <- aa_class_for_display(
          display_char = display_char,
          ref_char = ref_chars[pos],
          raw_char = raw_char,
          after_stop = after_stop[pos],
          is_ref = is_ref
        )
        
        span(
          class = css_class,
          title = paste0(
            "Sequence: ", seq_names[i],
            "\nAlignment column: ", pos,
            "\nReference AA: ", ref_chars[pos],
            "\nDisplayed AA: ", display_char,
            "\nRaw AA: ", raw_char
          ),
          display_char
        )
      })
      
      row_class <- if (is_ref) {
        "bioedit-row bioedit-ref-row"
      } else {
        "bioedit-row"
      }
      
      label_class <- if (is_ref) {
        "bioedit-label bioedit-ref-label"
      } else {
        "bioedit-label"
      }
      
      rows[[length(rows) + 1]] <- div(
        class = row_class,
        span(
          class = label_class,
          title = seq_names[i],
          seq_names[i]
        ),
        span(
          class = "bioedit-seq",
          aa_spans
        )
      )
    }
    
    tagList(rows)
  })
  
  output$download_aln_fasta <- downloadHandler(
    filename = function() {
      paste0(
        gene_data$current_gene,
        "_variant_protein_clustalomega_aln_",
        Sys.Date(),
        ".fasta"
      )
    },
    content = function(file) {
      req(gene_data$msa_result)
      req(gene_data$msa_result$alignment_set)
      
      Biostrings::writeXStringSet(
        gene_data$msa_result$alignment_set,
        filepath = file
      )
    }
  )
  
  output$download_genotypes <- downloadHandler(
    filename = function() {
      paste0(
        gene_data$current_gene,
        "_genotypes_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      gt <- filtered_genotypes()
      
      if (is.null(gt)) {
        write.csv(
          data.frame(Message = "No genotypes available."),
          file,
          row.names = FALSE
        )
      } else {
        write.csv(gt, file, row.names = FALSE)
      }
    }
  )
  
  output$download_vcf <- downloadHandler(
    filename = function() {
      paste0(
        gene_data$current_gene,
        "_filtered_",
        Sys.Date(),
        ".vcf"
      )
    },
    content = function(file) {
      req(gene_data$vcf)
      
      variants <- filtered_variants()
      
      if (nrow(variants) == 0) {
        stop("No variants available to write.")
      }
      
      vcf_pos <- vcfR::getFIX(gene_data$vcf)[, "POS"]
      
      variant_indices <- which(as.numeric(vcf_pos) %in% variants$POS)
      
      if (length(variant_indices) == 0) {
        stop("No matching variants found in VCF.")
      }
      
      vcf_filtered <- gene_data$vcf[variant_indices, ]
      
      vcfR::write.vcf(vcf_filtered, file = file)
    }
  )
}

# ===== Run the app =====
shinyApp(ui = ui, server = server)