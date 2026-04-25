# Genotype calling based on coverage breadth
# Hedhly A. et al. (2025)
# Input: bedtools coverage output (.tsv)
# Output: genotype table (CSV)
#
# Note: "V31" is the reference-version label used in Hedhly A. et al. (2025).
# If using another reference, update input paths, file names, and output labels
# accordingly.

############################################################
# User-defined paths and parameters
############################################################

data_dir <- "path/to/coverage/files"
output_dir <- "path/to/output"

# Minimum Breadth of coverage required for allele calling
min_breadth <- 0.50

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# Collect coverage files
############################################################

coverage_pattern <- "_F0x4q20eNM3_BedT_coverage_hist\\.tsv$"
coverage_files <- list.files(
  path = data_dir,
  pattern = coverage_pattern,
  full.names = TRUE
)

if (length(coverage_files) == 0) {
  stop("No coverage files found in: ", data_dir)
}

############################################################
# Functions
############################################################

extract_gene_id <- function(attributes) {
  gene <- sub(".*ID=([^;]+).*", "\\1", attributes)
  missing_id <- gene == attributes

  gene[missing_id] <- sub(
    ".*Name=([^;]+).*",
    "\\1",
    attributes[missing_id]
  )

  gene
}

read_coverage_file <- function(file) {
  df_coverage <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    na.strings = "NA",
    stringsAsFactors = FALSE
  )

  # Extract gene identifier from the GFF attribute field
  df_coverage$Gene <- extract_gene_id(as.character(df_coverage$V10))

  # Bedtools coverage fraction / breadth
  df_coverage$Breadth <- as.numeric(df_coverage$V14)

  df_coverage
}

call_genotype <- function(file, min_breadth = 0.50) {
  df_coverage <- read_coverage_file(file)

  # Filter by Breadth of coverage
  df_called <- df_coverage[df_coverage$Breadth >= min_breadth, ]
  df_called$Breadth <- round(df_called$Breadth, digits = 2)

  # Extract accession ID by removing pipeline-specific suffix
  accession <- sub(coverage_pattern, "", basename(file))

  if (nrow(df_called) == 0) {
    genotype <- NA_character_
    breadth <- NA_character_
  } else {
    genotype <- paste(df_called$Gene, collapse = "|")
    breadth <- paste(df_called$Breadth, collapse = "|")
  }

  data.frame(
    Accession = accession,
    Genotype = genotype,
    Breadth = breadth,
    stringsAsFactors = FALSE
  )
}

############################################################
# Run genotype calling
############################################################

genotypes_table <- do.call(
  rbind,
  lapply(coverage_files, call_genotype, min_breadth = min_breadth)
)

############################################################
# Save output
############################################################

write.csv(
  genotypes_table,
  file = file.path(output_dir, "Genotypes_table.csv"),
  row.names = FALSE
)
