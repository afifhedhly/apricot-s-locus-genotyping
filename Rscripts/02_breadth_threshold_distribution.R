# Breadth histogram with a manually split y-axis
# Hedhly A. et al. (2025)
# Input: bedtools coverage output (.tsv)
# Output: histogram plot (PNG/SVG) and plot data (CSV)
#
# Note: "V31" is the reference-version label used in this study.
# If using another reference, update input paths, file names, and output labels
# accordingly.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(svglite)
})

############################################################
# User-defined paths and parameters
############################################################

data_dir <- "path/to/coverage/files"
output_dir <- "path/to/output"

# Breadth threshold shown as a vertical dashed line
threshold <- 0.50

# Histogram and split-axis settings
bins <- 60

# Omitted interval in the y-axis. Adjust these two values to match
# the histogram counts in your dataset, or set them so high that the
# script falls back to a normal histogram without a y-axis break.
# The split is made by stacking two y-axis ranges, not by using a
# continuous axis-break package.
y_break <- c(1000, 12000)

output_prefix <- "Breadth_distribution"

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

read_breadth_values <- function(file) {
  df_coverage <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    na.strings = "NA",
    stringsAsFactors = FALSE
  )

  # V14 is the bedtools coverage fraction / breadth
  as.numeric(df_coverage$V14)
}

make_split_breadth_histogram <- function(df_breadth,
                                         threshold = 0.50,
                                         bins = 60,
                                         y_break = c(1000, 12000)) {
  hist_counts <- hist(df_breadth$Breadth, breaks = bins, plot = FALSE)$counts
  top_limit <- max(hist_counts, na.rm = TRUE) * 1.05

  base_plot <- ggplot(df_breadth, aes(x = Breadth)) +
    geom_histogram(
      bins = bins,
      color = "black",
      fill = "grey75",
      linewidth = 0.2
    ) +
    geom_vline(
      xintercept = threshold,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    labs(x = "Breadth of coverage", y = "Frequency") +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  if (top_limit <= y_break[2]) {
    return(
      base_plot +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, top_limit)) +
        labs(title = "Distribution of Breadth values") +
        theme(plot.title = element_text(face = "bold"))
    )
  }

  top_plot <- base_plot +
    coord_cartesian(xlim = c(0, 1), ylim = c(y_break[2], top_limit)) +
    labs(title = "Distribution of Breadth values") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )

  bottom_plot <- base_plot +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, y_break[1])) +
    theme(
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    )

  grid.arrange(top_plot, bottom_plot, ncol = 1, heights = c(1, 2))
}

############################################################
# Collect breadth values
############################################################

all_breadth <- unlist(lapply(coverage_files, read_breadth_values))

df_breadth <- data.frame(
  Breadth = all_breadth[!is.na(all_breadth)],
  stringsAsFactors = FALSE
)

write.csv(
  df_breadth,
  file = file.path(output_dir, paste0(output_prefix, "_plot_data.csv")),
  row.names = FALSE
)

############################################################
# Plot and save
############################################################

histogram_grob <- make_split_breadth_histogram(
  df_breadth = df_breadth,
  threshold = threshold,
  bins = bins,
  y_break = y_break
)

ggsave(
  file.path(output_dir, paste0(output_prefix, "_y_axis_break.png")),
  plot = histogram_grob,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(output_dir, paste0(output_prefix, "_y_axis_break.svg")),
  plot = histogram_grob,
  width = 7,
  height = 5,
  device = svglite
)
