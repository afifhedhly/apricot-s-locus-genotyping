# Genotype calling based on coverage breadth
# Lora J. et al. (2025)
# Input: bedtools coverage output (.tsv)
# Output: genotype table (CSV)
# Define directories
data_dir <- "path/to/coverage/files"
rDir <- "path/to/output"# Get a list of all tsv files generated with “bedtools coverage” command (one per accession),coverage_files <- list.files(path = data_dir, pattern = "_F0x4q20eNM3_BedT_coverage_hist\\.tsv$")genotypes_table = data.frame()

for(value in coverage_files) {

  file <- file.path(data_dir, value)

  df_coverage <- as.data.frame(read.table(file, header=F, sep='\t',
                                          na.strings="NA", stringsAsFactors=FALSE))

  # Extract gene names
  df_coverage$V10 <- sapply(strsplit(as.character(df_coverage$V10), '='), `[`, 2)

  # Filter by coverage
  df_coverage <- df_coverage[df_coverage$V14 >= 0.50, ]

  # Round coverage
  df_coverage$V14 <- round(df_coverage$V14, digits = 2)

  # Collapse results
  if (nrow(df_coverage) == 0) {
    genotype <- NA
    breadth <- NA
  } else {
    genotype <- paste(df_coverage$V10, collapse = "|")
    breadth <- paste(df_coverage$V14, collapse = "|")
  }

  # Extract accession ID by removing pipeline-specific suffix
  accession <- sub("_F0x4q20eNM3_BedT_coverage_hist.*", "", value)

  cols <- c(accession, genotype, breadth)

  genotypes_table = rbind(genotypes_table, cols)
}

# Assign column names ONCE
colnames(genotypes_table) <- c("Accession", "Genotype", "Breadth")

# Write output
write.csv(genotypes_table,
          file=file.path(rDir, "Genotypes_table.csv"),
          row.names = FALSE)