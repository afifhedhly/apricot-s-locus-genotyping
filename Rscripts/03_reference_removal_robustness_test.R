# Reference-removal robustness test
# Hedhly A. et al. (2025)
# Input: bedtools coverage output from full and modified references
# Output: robustness plot (PNG/PDF/SVG) and summary tables (CSV)
#
# Note: "V31" is the reference-version label used in this study.
# If using another reference, update input paths, file names, reference BED
# files, and output labels accordingly.

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

############################################################
# User-defined paths and parameters
############################################################

# Folder containing coverage files generated against the full reference
full_cov_dir <- "path/to/full_reference/coverage/files"

# Folder containing coverage files generated after remapping to modified
# references lacking the target allele(s)
remap_cov_dir <- "path/to/modified_reference/coverage/files"

# BED file for the full reference. This is used only to order genes in the plot.
full_ref_bed <- "path/to/full_reference_genes.bed"

output_dir <- "path/to/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Breadth threshold shown as a horizontal dashed line
breadth_threshold <- 0.50

# Genes with Breadth >= this value are retained for plotting, in addition
# to target removed genes and the top genes per reference.
plot_min_breadth <- 0.05

# Number of top Breadth genes retained per case/reference
top_n_genes <- 10

############################################################
# Define remapping cases
############################################################

# RemapCaseName must match the label used in modified-reference file names:
# <Accession>_case_<RemapCaseName>_F0x4q20eNM3_BedT_coverage_hist.tsv
cases <- data.frame(
  Case = c("S13 removal", "S6 removal", "S4 + S23 removal"),
  Accession = c("ERR5013362", "ERR4762650", "ERR4762463"),
  RemapCaseName = c("S13", "S6", "S4S23"),
  stringsAsFactors = FALSE
)

# Target genes are the BED Name= values missing from each modified reference.
target_genes <- list(
  "S13 removal" = c("SFB5", "S13"),
  "S6 removal" = c("SFB6", "S6"),
  "S4 + S23 removal" = c("SFB4", "S4", "SFB23", "S23")
)

############################################################
# Helper functions
############################################################

extract_gene_id <- function(attributes) {
  gene <- sub(".*Name=([^;]+).*", "\\1", attributes)
  missing_name <- gene == attributes

  gene[missing_name] <- sub(
    ".*ID=([^;]+).*",
    "\\1",
    attributes[missing_name]
  )

  gene <- sub("\r$", "", gene)
  gene
}

read_bedtools_coverage <- function(file, case_name, accession, reference_type) {
  if (!file.exists(file)) {
    stop("Missing coverage file: ", file)
  }

  tab <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    na.strings = "NA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(tab) < 14) {
    stop("Expected at least 14 columns in: ", file)
  }

  breadth <- as.numeric(tab[[14]])

  if (any(is.na(breadth))) {
    warning("NA breadth values found in: ", file)
  }

  data.frame(
    Case = case_name,
    Accession = accession,
    Reference = reference_type,
    Gene = extract_gene_id(as.character(tab[[10]])),
    Breadth = breadth,
    stringsAsFactors = FALSE
  )
}

top_genes_by_reference <- function(dat, n = 10) {
  by_reference <- split(dat, dat$Reference)

  unique(unlist(lapply(by_reference, function(ref_dat) {
    ref_dat <- ref_dat[order(ref_dat$Breadth, decreasing = TRUE), , drop = FALSE]
    head(ref_dat$Gene, n)
  }), use.names = FALSE))
}

read_bed_gene_order <- function(file) {
  if (!file.exists(file)) {
    stop("Missing full-reference BED file: ", file)
  }

  bed <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  unique(extract_gene_id(as.character(bed[[10]])))
}

fallback_gene_order <- function(genes) {
  genes <- unique(as.character(genes))

  sort_key <- lapply(genes, function(gene) {
    is_sn <- grepl("^Sn[0-9]+$", gene)

    if (is_sn) {
      return(data.frame(
        Gene = gene,
        Group = 2L,
        Allele = as.integer(sub("^Sn", "", gene)),
        Variant = 0L,
        Type = 0L,
        Label = gene,
        stringsAsFactors = FALSE
      ))
    }

    type <- if (grepl("^SFB", gene)) {
      1L
    } else if (grepl("^S", gene)) {
      0L
    } else {
      2L
    }

    suffix <- sub("^SFB", "", gene)
    suffix <- sub("^S", "", suffix)
    allele <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", suffix)))

    if (is.na(allele)) {
      allele <- .Machine$integer.max
    }

    h_variant <- grepl("H[0-9]+", suffix)
    variant <- if (h_variant) {
      h_number <- suppressWarnings(as.integer(sub(".*H([0-9]+).*", "\\1", suffix)))
      100L + ifelse(is.na(h_number), 0L, h_number)
    } else if (grepl("[A-Za-z]", suffix)) {
      10L
    } else {
      0L
    }

    data.frame(
      Gene = gene,
      Group = 1L,
      Allele = allele,
      Variant = variant,
      Type = type,
      Label = gene,
      stringsAsFactors = FALSE
    )
  })

  sort_key <- do.call(rbind, sort_key)
  sort_key$Gene[order(
    sort_key$Group,
    sort_key$Allele,
    sort_key$Variant,
    sort_key$Type,
    sort_key$Label
  )]
}

gene_order_levels <- function(genes, reference_order) {
  genes <- unique(as.character(genes))

  ordered_genes <- reference_order[reference_order %in% genes]
  missing_from_reference <- setdiff(genes, ordered_genes)

  if (length(missing_from_reference) > 0) {
    ordered_genes <- c(ordered_genes, fallback_gene_order(missing_from_reference))
  }

  ordered_genes
}

############################################################
# Read full-reference and remapped coverage files
############################################################

all_cov <- data.frame()

for (i in seq_len(nrow(cases))) {
  case_name <- cases$Case[i]
  acc <- cases$Accession[i]
  remap_case <- cases$RemapCaseName[i]

  full_file <- file.path(
    full_cov_dir,
    paste0(acc, "_F0x4q20eNM3_BedT_coverage_hist.tsv")
  )

  remap_file <- file.path(
    remap_cov_dir,
    paste0(acc, "_case_", remap_case, "_F0x4q20eNM3_BedT_coverage_hist.tsv")
  )

  all_cov <- rbind(
    all_cov,
    read_bedtools_coverage(
      file = full_file,
      case_name = case_name,
      accession = acc,
      reference_type = "Full reference"
    ),
    read_bedtools_coverage(
      file = remap_file,
      case_name = case_name,
      accession = acc,
      reference_type = "Reference lacking target allele(s)"
    )
  )
}

# Collapse duplicated gene labels per case/reference so plotting does not
# overdraw bars when repeated labels exist in the reference.
all_cov <- aggregate(
  Breadth ~ Case + Accession + Reference + Gene,
  data = all_cov,
  FUN = max
)

############################################################
# Select informative genes for plotting
############################################################

selected_genes <- c()

for (case_name in unique(all_cov$Case)) {
  case_data <- all_cov[all_cov$Case == case_name, , drop = FALSE]

  target <- target_genes[[case_name]]
  nonzero_genes <- case_data$Gene[case_data$Breadth >= plot_min_breadth]
  top_genes <- top_genes_by_reference(case_data, n = top_n_genes)

  selected_genes <- c(selected_genes, target, nonzero_genes, top_genes)
}

selected_genes <- unique(selected_genes)
plot_data <- all_cov[all_cov$Gene %in% selected_genes, , drop = FALSE]

############################################################
# Add missing target genes in modified references
############################################################

# When an allele is removed from the modified reference, it will not exist
# in the remap BED file. For visualization, add it with Breadth = 0.
missing_rows <- data.frame()

for (case_name in names(target_genes)) {
  acc <- cases$Accession[cases$Case == case_name]

  for (gene in target_genes[[case_name]]) {
    exists_in_remap <- any(
      plot_data$Case == case_name &
        plot_data$Reference == "Reference lacking target allele(s)" &
        plot_data$Gene == gene
    )

    if (!exists_in_remap) {
      missing_rows <- rbind(
        missing_rows,
        data.frame(
          Case = case_name,
          Accession = acc,
          Reference = "Reference lacking target allele(s)",
          Gene = gene,
          Breadth = 0,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

plot_data <- rbind(plot_data, missing_rows)

plot_data$TargetRemoved <- ifelse(
  mapply(
    function(case, gene) gene %in% target_genes[[case]],
    plot_data$Case,
    plot_data$Gene
  ),
  "Target allele",
  "Other allele"
)

plot_data$Case <- factor(plot_data$Case, levels = cases$Case)
plot_data$Reference <- factor(
  plot_data$Reference,
  levels = c("Full reference", "Reference lacking target allele(s)")
)

full_ref_gene_order <- read_bed_gene_order(full_ref_bed)
plot_data$Gene <- factor(
  plot_data$Gene,
  levels = gene_order_levels(plot_data$Gene, full_ref_gene_order)
)
plot_data <- plot_data[
  order(plot_data$Case, plot_data$Gene, plot_data$Reference),
  ,
  drop = FALSE
]

############################################################
# Save summary tables
############################################################

write.csv(
  plot_data,
  file = file.path(output_dir, "remapping_breadth_plot_data.csv"),
  row.names = FALSE
)

nonzero_summary <- all_cov[all_cov$Breadth > 0, , drop = FALSE]
nonzero_summary <- nonzero_summary[
  order(nonzero_summary$Case, nonzero_summary$Reference, -nonzero_summary$Breadth),
  ,
  drop = FALSE
]

write.csv(
  nonzero_summary,
  file = file.path(output_dir, "remapping_all_nonzero_breadth.csv"),
  row.names = FALSE
)

############################################################
# Plot and save
############################################################

p <- ggplot(
  plot_data,
  aes(x = Gene, y = Breadth, fill = Reference)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_hline(
    yintercept = breadth_threshold,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(~ Case, scales = "free_x", ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Effect of removing target alleles from the synthetic reference",
    x = "Reference allele / gene",
    y = "Breadth of coverage",
    fill = "Reference"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p)

ggsave(
  filename = file.path(output_dir, "remapping_breadth_robustness.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "remapping_breadth_robustness.pdf"),
  plot = p,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(output_dir, "remapping_breadth_robustness.svg"),
  plot = p,
  width = 10,
  height = 8,
  device = svglite
)
