#!/usr/bin/env Rscript

###############################################################################
# run_CNAscreen.R
#
# Purpose:
#   Run a QDNAseq-based genome integrity screen from BAM files and produce
#   genome-wide copy-number plots suitable for cultured-cell QC
#   (e.g. iPSCs, fibroblasts, ESCs).
#
# Main outputs per sample:
#   - <sample>.plot.png
#   - <sample>.qdnaseq.rds
#   - <sample>.bins.tsv.gz
#   - <sample>.segments.tsv
#
# Run-level outputs:
#   - sample_qc_summary.tsv
#   - SESSION_INFO.txt
#   - CNAscreen_summary.pdf
# Notes:
#   - Copy number is plotted as log2 ratio:
#       ~0   = copy-number neutral / diploid baseline
#       >0   = relative gain
#       <0   = relative loss
#   - Intended for broad genome integrity QC, not focal CNV discovery.
#
# Runs from any shell (a terminal, an RStudio Terminal tab, an HPC job
# script) - it has no RStudio dependency:
#   Rscript run_CNAscreen.R
#     --bam-dir bam
#     --out-dir CNAscreen_output
#     --bin-size 100
#     --genome hg38
#     --sex-map sample_sex.csv
#     --map-cutoff 80
#     --bases-cutoff 99.9
#     --seed 170826
#
###############################################################################

if (!startsWith(as.character(getRversion()), "4.5.")) {
  stop("CNAscreen requires R 4.5.x, but this is R ", getRversion(), ". ",
       "See README > Installation for how to install R 4.5.3 with rig.",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(optparse)
  library(tools)
  library(grid)
  library(BiocParallel)
})

SCRIPT_VERSION <- "1.1.0"

# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------

option_list <- list(
  make_option(
    "--bam-dir",
    dest = "bam_dir",
    type = "character",
    default = "bam",
    help = "Directory containing BAM (.bam) files and their indexes (.bam.bai) [default: %default]"
  ),
  make_option(
    "--out-dir",
    dest = "out_dir",
    type = "character",
    default = "CNAscreen_output",
    help = "Output directory [default: %default]"
  ),
  make_option(
    "--bin-size",
    dest = "bin_size",
    type = "integer",
    default = 100,
    help = "Bin size in kb for QDNAseq [default: %default]"
  ),
  make_option(
    "--genome",
    dest = "genome",
    type = "character",
    default = "hg38",
    help = "Genome build used for bin annotations, i.e. hg38 [default: %default]"
  ),
  make_option(
    "--sex-map",
    dest = "sex_map",
    type = "character",
    default = NULL,
    help = "REQUIRED. CSV with columns sample,sex where sex is XX or XY. Every sample in --bam-dir must have an entry."
  ),
  make_option(
    "--no-sex-ploidy-correction",
    dest = "no_sex_ploidy_correction",
    action = "store_true",
    default = FALSE,
    help = paste(
      "Disable the XY sex-chromosome ploidy correction (rescales chrX/chrY",
      "copy number by 2x for XY samples so they sit at the same baseline",
      "as autosomes). The correction is applied by default; set this flag",
      "to see the raw, uncorrected QDNAseq output instead [default: %default]"
    )
  ),
  make_option(
    "--map-cutoff",
    dest = "map_cutoff",
    type = "double",
    default = 80,
    help = "Minimum mappability percentage for QC flagging [default: %default]"
  ),
  make_option(
    "--bases-cutoff",
    dest = "bases_cutoff",
    type = "double",
    default = 99.9,
    help = "Minimum characterized bases percentage for QC flagging [default: %default]"
  ),
  make_option(
    "--exclude-regions",
    dest = "exclude_regions",
    type = "character",
    default = NULL,
    help = paste(
      "Optional BED file (chrom, start, end; tab-delimited, 0-based",
      "half-open per BED convention; extra columns ignored; 'chr' prefix",
      "on chromosome names is stripped automatically) of additional",
      "genomic regions to mask, beyond QDNAseq's built-in blacklist,",
      "mappability, and bases filters. Use this for recurrent",
      "artifact bins identified empirically (e.g. via cross-sample",
      "comparison) or known problem regions such as segmental",
      "duplications, e.g. exported directly from the UCSC Table Browser",
      "or ENCODE blacklist. Bins overlapping any listed region are",
      "flagged with status 'Excluded region' and excluded from the",
      "genome-wide plot and QC summary the same way as other failing",
      "bins [default: none]"
    )
  ),
  make_option(
    "--min-mapq",
    dest = "min_mapq",
    type = "integer",
    default = 37,
    help = "Minimum BAM mapping quality passed to binReadCounts() [default: %default]"
  ),
  make_option(
    "--segment-alpha",
    dest = "segment_alpha",
    type = "double",
    default = 1e-10,
    help = paste(
      "DNAcopy significance level for accepting a segment change-point in",
      "segmentBins(). Smaller = fewer, more conservative segments. The",
      "default is QDNAseq's own. Looser values can over-segment low-noise samples into many",
      "small segments; stricter values reduce sensitivity to small or",
      "subclonal events [default: %default]"
    )
  ),
  make_option(
    "--segment-undo-sd",
    dest = "segment_undo_sd",
    type = "double",
    default = 1,
    help = paste(
      "DNAcopy undo.SD used by segmentBins() with undo.splits = 'sdundo':",
      "adjacent segments whose means differ by less than this many SDs of",
      "the data are merged. Because it is in SD units it scales with each",
      "sample's noise. The default is QDNAseq's own [default: %default]"
    )
  ),
  make_option(
    "--segment-tol",
    dest = "segment_tol",
    type = "double",
    default = 0.01,
    help = "Tolerance for merging adjacent equal segment values in plotting [default: %default]"
  ),
  make_option(
    "--seed",
    dest = "seed",
    type = "integer",
    default = 170826,
    help = "Random seed for reproducible DNAcopy segmentation [default: %default]"
  ),
  make_option(
    "--plot-filtered-bins",
    dest = "plot_filtered_bins",
    type = "logical",
    default = TRUE,
    help = "Show filtered (low-quality) bins as light grey points: TRUE or FALSE [default: %default]"
  ),
  make_option(
    "--width-mm",
    dest = "width_mm",
    type = "double",
    default = 95,
    help = "Plot width in mm [default: %default]"
  ),
  make_option(
    "--height-mm",
    dest = "height_mm",
    type = "double",
    default = 60,
    help = "Plot height in mm [default: %default]"
  ),
  make_option(
    "--ymin",
    dest = "ymin",
    type = "double",
    default = -2.2,
    help = "Plot lower y-limit [default: %default]"
  ),
  make_option(
    "--ymax",
    dest = "ymax",
    type = "double",
    default = 2.2,
    help = "Plot upper y-limit [default: %default]"
  ),
  make_option(
    "--verbose",
    dest = "verbose",
    type = "logical",
    default = TRUE,
    help = "Print progress messages: TRUE or FALSE [default: %default]"
  ),
  make_option(
    "--threads",
    dest = "threads",
    type = "integer",
    default = 1,
    help = paste(
      "Number of samples to process in parallel [default: %default,",
      "serial]. Size this against your available RAM and BAM file sizes,",
      "not sample count: each worker holds one sample's read data in",
      "memory at once, so N workers can use roughly N times the memory",
      "of a single serial run. A handful of large/high-coverage BAMs",
      "run in parallel can exhaust a laptop's RAM faster than running",
      "them serially would. Uses BiocParallel::MulticoreParam (forked",
      "workers) on Mac/Linux. Windows does not support fork-based",
      "workers, so on Windows this flag is ignored above 1 and the run",
      "falls back to serial execution with a warning, rather than use",
      "an unverified separate-process backend (BiocParallel::SnowParam)",
      "that could silently produce incorrect results; run under WSL for",
      "parallel throughput on Windows. Capped automatically at the",
      "number of samples and at (detected CPU cores - 1)."
    )
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list
)

opt <- parse_args(parser)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

msg <- function(...) {
  if (isTRUE(opt$verbose)) {
    message(...)
  }
}

stop_if_not_dir <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !dir.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
}

validate_writable_outdir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  testfile <- file.path(path, paste0(".write_test_", Sys.getpid()))
  ok <- tryCatch(
    {
      file.create(testfile)
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  if (!isTRUE(ok)) {
    stop("Output directory is not writable: ", path, call. = FALSE)
  }
  unlink(testfile)
}

read_sex_map <- function(path) {
  if (!file.exists(path)) {
    stop("Sex map file not found: ", path, call. = FALSE)
  }
  
  x <- read.csv(path, stringsAsFactors = FALSE)
  
  if (!all(c("sample", "sex") %in% names(x))) {
    stop("Sex map must contain columns: sample, sex", call. = FALSE)
  }
  
  x$sample <- trimws(gsub("^\ufeff", "", as.character(x$sample)))
  x$sex <- toupper(trimws(as.character(x$sex)))
  
  bad <- setdiff(unique(x$sex), c("XX", "XY"))
  if (length(bad) > 0) {
    stop("Invalid sex values in sex map: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  
  if (anyDuplicated(x$sample)) {
    dup <- unique(x$sample[duplicated(x$sample)])
    stop("Duplicate sample names in sex map: ", paste(dup, collapse = ", "), call. = FALSE)
  }
  
  x
}

safe_log2 <- function(x, min_positive = 1e-3) {
  ifelse(is.finite(x) & x >= min_positive, log2(x), NA_real_)
}

read_exclude_regions <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("--exclude-regions file not found: ", path, call. = FALSE)
  }

  raw <- tryCatch(
    read.table(
      path,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      comment.char = "#"
    ),
    error = function(e) {
      stop(
        "Could not parse --exclude-regions file as tab-delimited BED '",
        path, "': ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (ncol(raw) < 3) {
    stop(
      "--exclude-regions file must be tab-delimited with at least 3 columns ",
      "(chrom, start, end): ", path,
      call. = FALSE
    )
  }

  regions <- tibble(
    chrom = gsub("^chr", "", as.character(raw[[1]]), ignore.case = TRUE),
    # BED is 0-based, half-open; convert to 1-based inclusive to match
    # QDNAseq bin coordinates (fData()$start / $end).
    start = as.numeric(raw[[2]]) + 1,
    end   = as.numeric(raw[[3]])
  )

  if (anyNA(regions$start) || anyNA(regions$end)) {
    stop(
      "--exclude-regions file has non-numeric start/end values: ", path,
      call. = FALSE
    )
  }
  if (any(regions$end < regions$start)) {
    stop(
      "--exclude-regions file contains rows where end < start: ", path,
      call. = FALSE
    )
  }

  regions
}

flag_excluded_bins <- function(df, regions) {
  if (is.null(regions) || nrow(regions) == 0) {
    df$excluded_region <- FALSE
    return(df)
  }

  bin_gr <- GenomicRanges::GRanges(
    seqnames = as.character(df$chromosome),
    ranges = IRanges::IRanges(start = df$start, end = df$end)
  )
  region_gr <- GenomicRanges::GRanges(
    seqnames = regions$chrom,
    ranges = IRanges::IRanges(start = regions$start, end = regions$end)
  )

  hits <- GenomicRanges::findOverlaps(bin_gr, region_gr)

  df$excluded_region <- FALSE
  df$excluded_region[S4Vectors::queryHits(hits)] <- TRUE
  df
}

chr_levels <- function(chr_vec, include_Y = TRUE) {
  lev <- c(as.character(1:22), "X", if (isTRUE(include_Y)) "Y" else NULL)
  factor(as.character(chr_vec), levels = lev[lev %in% unique(as.character(chr_vec))])
}

mat_to_long <- function(mat, value_name) {
  df <- cbind(bin = seq_len(nrow(mat)), as.data.frame(mat, check.names = FALSE))
  pivot_longer(df, -bin, names_to = "sample", values_to = value_name)
}

apply_sex_ploidy_correction <- function(df, is_xy = FALSE) {
  if (!isTRUE(is_xy)) {
    return(df)
  }
  df |>
    mutate(
      CN_raw = ifelse(chromosome %in% c("X", "Y"), CN_raw * 2, CN_raw),
      SEG_raw = ifelse(chromosome %in% c("X", "Y"), SEG_raw * 2, SEG_raw)
    )
}

extract_plot_df <- function(cn, include_Y = TRUE, is_xy = include_Y) {
  fd <- as.data.frame(Biobase::fData(cn))
  fd$bin <- seq_len(nrow(fd))
  
  bases_pct <- as.numeric(fd$bases)
  map_raw <- as.numeric(fd$mappability)
  
  mx <- suppressWarnings(max(map_raw, na.rm = TRUE))
  mappability_pct <- if (is.finite(mx) && mx <= 1.1) 100 * map_raw else map_raw
  mappability_pct <- pmin(pmax(mappability_pct, 0), 100)
  
  bins <- tibble(
    bin = fd$bin,
    chromosome = chr_levels(fd$chromosome, include_Y = include_Y),
    start = as.numeric(fd$start),
    end = as.numeric(fd$end),
    bases_pct = bases_pct,
    mappability_pct = mappability_pct
  ) |>
    filter(!is.na(chromosome))
  
  assays <- Biobase::assayDataElementNames(cn)
  
  cn_long <- mat_to_long(Biobase::assayDataElement(cn, "copynumber"), "CN_raw")
  
  seg_long <- if ("segmented" %in% assays) {
    mat_to_long(Biobase::assayDataElement(cn, "segmented"), "SEG_raw")
  } else {
    tibble(bin = integer(), sample = character(), SEG_raw = numeric())
  }
  
  out <- bins |>
    left_join(cn_long, by = "bin") |>
    left_join(seg_long, by = c("bin", "sample")) |>
    arrange(sample, chromosome, start) |>
    as_tibble() |>
    apply_sex_ploidy_correction(is_xy = is_xy)

  out |>
    mutate(
      CN_log2  = safe_log2(CN_raw),
      SEG_log2 = safe_log2(SEG_raw)
    )
}

add_qc_flags <- function(df, bases_cut, map_cut) {
  if (!"excluded_region" %in% names(df)) {
    df$excluded_region <- FALSE
  }
  df |>
    mutate(
      low_bases = is.finite(bases_pct) & bases_pct < bases_cut,
      low_map = is.finite(mappability_pct) & mappability_pct < map_cut,
      status = case_when(
        excluded_region & (low_bases | low_map) ~ "Excluded + QC fail",
        excluded_region ~ "Excluded region",
        low_bases & low_map ~ "Both low",
        low_bases ~ "Bases < cutoff",
        low_map ~ "Low mappability",
        TRUE ~ "OK"
      )
    )
}

get_chr_offsets <- function(cn, include_Y = TRUE) {
  fa <- as.data.frame(Biobase::fData(cn))
  fa$chromosome <- chr_levels(fa$chromosome, include_Y = include_Y)
  fa <- subset(fa, !is.na(chromosome))
  fa$end <- as.numeric(fa$end)
  
  fa |>
    group_by(chromosome) |>
    summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") |>
    arrange(chromosome) |>
    mutate(
      chr_len = as.numeric(chr_len),
      offset = lag(cumsum(chr_len), default = 0),
      xmin_mb = offset / 1e6,
      xmax_mb = (offset + chr_len) / 1e6,
      center_mb = (offset + chr_len / 2) / 1e6,
      band = (row_number() %% 2) == 1
    )
}

add_pos_cum <- function(df_in, offs) {
  df_in |>
    mutate(
      chromosome = factor(as.character(chromosome), levels = levels(offs$chromosome)),
      start = as.numeric(start),
      end = as.numeric(end)
    ) |>
    left_join(select(offs, chromosome, offset), by = "chromosome") |>
    mutate(
      pos_cum = (start + offset) / 1e6,
      pos_cum_end = (end + offset) / 1e6
    ) |>
    arrange(chromosome, start)
}

compute_segment_runs <- function(seg_df, tol = 1e-3) {
  seg_df |>
    filter(
      is.finite(SEG_log2),
      is.finite(pos_cum),
      is.finite(pos_cum_end)
    ) |>
    arrange(chromosome, start) |>
    group_by(chromosome) |>
    mutate(
      new_run = row_number() == 1 | abs(SEG_log2 - lag(SEG_log2)) > tol
    ) |>
    mutate(run_id = cumsum(new_run)) |>
    group_by(chromosome, run_id) |>
    summarise(
      start   = first(start),
      end     = last(end),
      n_bins  = n(),
      seg_y   = first(SEG_log2),
      x_minMb = first(pos_cum),
      x_maxMb = last(pos_cum_end),
      .groups = "drop"
    )
}

limited_chr_labels <- function(offs) {
  wanted <- c(as.character(1:10), "12", "14", "16", "18", "21", "X", "Y")
  labs <- as.character(offs$chromosome)
  labs[!(labs %in% wanted)] <- ""
  list(breaks = offs$center_mb, labels = labs)
}

plot_genomewide <- function(ok_df, filtered_df, seg_runs, offs, ymin = -2.2, ymax = 2.2,
                            plot_filtered_bins = FALSE) {
  total_len_mb <- max(offs$xmax_mb, na.rm = TRUE)
  lab <- limited_chr_labels(offs)
  chr_boundaries <- offs$xmax_mb[-nrow(offs)]
  
  p <- ggplot() +
    geom_rect(
      data = offs,
      aes(xmin = xmin_mb, xmax = xmax_mb, ymin = -Inf, ymax = Inf, fill = band),
      inherit.aes = FALSE,
      alpha = 0.05,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = chr_boundaries,
      linewidth = 0.25,
      colour = "grey70"
    ) +
    geom_point(
      data = ok_df,
      aes(x = pos_cum, y = CN_log2),
      colour = "black",
      alpha = 1,
      size = 0.40,
      shape = 16,
      stroke = 0
    )
  
  if (isTRUE(plot_filtered_bins) && nrow(filtered_df) > 0) {
    p <- p +
      geom_point(
        data = filtered_df,
        aes(x = pos_cum, y = CN_log2),
        colour = "grey70",
        alpha = 0.4,
        size = 0.35,
        shape = 16,
        stroke = 0
      )
  }
  
  p +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "seagreen3"
    ) +
    geom_segment(
      data = seg_runs,
      aes(x = x_minMb, xend = x_maxMb, y = seg_y, yend = seg_y),
      colour = "#6A0DAD",
      linewidth = 1.6,
      alpha = 0.98,
      lineend = "butt"
    ) +
    scale_fill_manual(values = c("TRUE" = "grey80", "FALSE" = "white")) +
    scale_y_continuous(
      "Relative copy number (log2 ratio)",
      breaks = seq(-2, 2, by = 0.5),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_x_continuous(
      "Chromosome",
      breaks = lab$breaks,
      labels = lab$labels,
      expand = expansion(mult = c(0, 0)),
      limits = c(0, total_len_mb)
    ) +
    coord_cartesian(ylim = c(ymin, ymax), clip = "on") +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks.length = unit(2, "pt"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 9, margin = margin(t = 4)),
      axis.title.y = element_text(size = 9, margin = margin(r = 6)),
      plot.margin = margin(3, 6, 2, 6)
    )
}

write_run_summary_pdf <- function(plots, path, per_page = 8,
                                  page_w_mm = 210, page_h_mm = 297) {
  if (!length(plots)) return(invisible(NULL))
  
  grDevices::pdf(path, width = page_w_mm / 25.4, height = page_h_mm / 25.4,
                 onefile = TRUE, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  
  pages <- split(seq_along(plots), ceiling(seq_along(plots) / per_page))
  
  for (pg in seq_along(pages)) {
    grid::grid.newpage()
    # Header row + a fixed number of panel rows, so panels are the same
    # size on every page, including a part-filled last page
    lay <- grid::grid.layout(
      nrow = per_page + 1, ncol = 1,
      heights = grid::unit(c(8, rep(1, per_page)), c("mm", rep("null", per_page)))
    )
    grid::pushViewport(grid::viewport(
      layout = lay,
      width  = grid::unit(page_w_mm - 16, "mm"),
      height = grid::unit(page_h_mm - 16, "mm")
    ))
    
    grid::grid.text(
      sprintf("CNAscreen run summary: relative copy number (log2 ratio)    page %d of %d",
              pg, length(pages)),
      vp = grid::viewport(layout.pos.row = 1),
      gp = grid::gpar(fontsize = 8)
    )
    
    for (k in seq_along(pages[[pg]])) {
      i <- pages[[pg]][k]
      p <- plots[[i]] +
        ggplot2::labs(title = names(plots)[i]) +
        ggplot2::theme(
          plot.title  = ggplot2::element_text(size = 7, face = "bold",
                                              margin = ggplot2::margin(b = 1)),
          axis.title  = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(size = 6),
          axis.text.y = ggplot2::element_text(size = 6)
        )
      print(p, vp = grid::viewport(layout.pos.row = k + 1, layout.pos.col = 1))
    }
    grid::popViewport()
  }
  invisible(path)
}

# Bin-to-bin noise: median absolute difference in log2 ratio between
# consecutive passing (status == "OK") bins on the same chromosome
# (MAPD-style; unscaled). Lower = cleaner. It is insensitive to genuine
# copy-number segments (a step change affects only one pair of bins), so it
# measures bin-level noise rather than biology, and is comparable across
# samples at the same --bin-size. "Consecutive" means consecutive among
# passing bins, so bins removed by filtering are bridged over.
compute_noise_mapd <- function(df) {
  d <- df |>
    filter(status == "OK", is.finite(CN_log2)) |>
    arrange(chromosome, start) |>
    group_by(chromosome) |>
    mutate(step = abs(CN_log2 - lag(CN_log2))) |>
    ungroup() |>
    pull(step)

  if (!any(is.finite(d))) {
    return(NA_real_)
  }
  median(d, na.rm = TRUE)
}

sample_qc_summary <- function(df, total_reads = NA_real_, seg_runs = NULL) {
  tibble(
    total_reads = as.numeric(total_reads),
    noise_mapd = compute_noise_mapd(df),
    n_segments = if (is.null(seg_runs)) NA_integer_ else nrow(seg_runs),
    n_bins = nrow(df),
    n_finite = sum(is.finite(df$CN_log2)),
    frac_negative = mean(df$CN_log2 < 0, na.rm = TRUE),
    median_log2 = median(df$CN_log2, na.rm = TRUE),
    median_bases_pct = median(df$bases_pct, na.rm = TRUE),
    median_mappability_pct = median(df$mappability_pct, na.rm = TRUE),
    frac_low_bases = mean(df$low_bases, na.rm = TRUE),
    frac_low_map = mean(df$low_map, na.rm = TRUE),
    frac_excluded_region = mean(df$excluded_region, na.rm = TRUE)
  )
}

write_session_info <- function(out_dir, seed) {
  lines <- c(
    paste("Script version:", SCRIPT_VERSION),
    paste("Date:", as.character(Sys.time())),
    paste("Command:", paste(commandArgs(), collapse = " ")),
    paste("Random seed:", seed),
    "",
    capture.output(sessionInfo())
  )
  writeLines(lines, con = file.path(out_dir, "SESSION_INFO.txt"))
}

validate_inputs <- function(bamfiles, bamnames, sex_map) {
  if (!length(bamfiles)) {
    stop("No BAM files found.", call. = FALSE)
  }
  
  if (anyDuplicated(bamnames)) {
    dup <- unique(bamnames[duplicated(bamnames)])
    stop(
      "Duplicate sample names inferred from BAM filenames: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }
  
  missing_bai <- bamfiles[!file.exists(paste0(bamfiles, ".bai"))]
  if (length(missing_bai)) {
    stop(
      "Missing BAM index (.bai) for: ",
      paste(basename(missing_bai), collapse = ", "),
      call. = FALSE
    )
  }
  
  not_found <- setdiff(sex_map$sample, bamnames)
  if (length(not_found)) {
    stop(
      "Sample(s) in sex map do not match any BAM-derived sample name: ",
      paste(not_found, collapse = ", "), ". ",
      "BAM-derived sample names are: ", paste(bamnames, collapse = ", "), ". ",
      call. = FALSE
    )
  }
  
  missing_from_map <- setdiff(bamnames, sex_map$sample)
  if (length(missing_from_map)) {
    stop(
      "Sample(s) found in BAM directory but missing from sex map: ",
      paste(missing_from_map, collapse = ", "), ". ",
      "Every BAM sample must have a sex assignment in the sex-map CSV.",
      call. = FALSE
    )
  }
}

export_table_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ---------------------------------------------------------------------------
# Input discovery and validation
# ---------------------------------------------------------------------------

stop_if_not_dir(opt$bam_dir, "BAM directory")
validate_writable_outdir(opt$out_dir)

if (is.null(opt$sex_map) || !nzchar(opt$sex_map)) {
  stop(
    "--sex-map is required. Provide a CSV with columns sample,sex (XX or XY) ",
    "covering every sample in --bam-dir.",
    call. = FALSE
  )
}

sex_map <- read_sex_map(opt$sex_map)

bam_dir_abs <- normalizePath(opt$bam_dir, winslash = "/", mustWork = TRUE)
out_dir_abs <- normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE)

bamfiles <- list.files(
  path = bam_dir_abs,
  pattern = "\\.bam$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)
bamfiles <- sort(bamfiles, method = "radix")

bamnames <- file_path_sans_ext(basename(bamfiles))

validate_inputs(bamfiles, bamnames, sex_map)

msg("Found ", length(bamfiles), " BAM file(s).")
msg("BAM directory: ", bam_dir_abs)
msg("Samples: ", paste(bamnames, collapse = ", "))
msg("Genome: ", opt$genome, " | Bin size: ", opt$bin_size, " kb")
msg("Output directory: ", out_dir_abs)

# ---------------------------------------------------------------------------
# Load bin annotations
# ---------------------------------------------------------------------------

msg("Loading bin annotations: genome=", opt$genome, ", bin-size=", opt$bin_size, " kb")
bins <- tryCatch(
  getBinAnnotations(binSize = opt$bin_size, genome = opt$genome),
  error = function(e) {
    stop(
      "Could not load QDNAseq bin annotations for genome='", opt$genome,
      "', bin-size=", opt$bin_size, " kb. ",
      "Open CNAscreen.Rproj and run renv::restore() in the Console.",
      call. = FALSE
    )
  }
)

# ---------------------------------------------------------------------------
# Per-sample processing
#
# process_one_sample() runs the full pipeline for one BAM: binning, filtering,
# correction, segmentation, plot and output files. It returns that sample's
# QC row and plot for the run-level summary. Samples are independent, so
# output is identical whether they run serially or in parallel.
# ---------------------------------------------------------------------------

exclude_regions <- read_exclude_regions(opt$exclude_regions)
if (!is.null(exclude_regions)) {
  msg(
    "Loaded ", nrow(exclude_regions), " exclude region(s) from ",
    opt$exclude_regions
  )
}

process_one_sample <- function(i, bamfiles, bamnames, sex_map, bins,
                                exclude_regions, opt) {
  # Loaded here so the function also works in a fresh worker session.
  suppressPackageStartupMessages({
    library(QDNAseq); library(Biobase); library(dplyr); library(tidyr)
    library(tibble); library(ggplot2)
  })

  bamfile <- bamfiles[i]
  sample_name <- bamnames[i]

  # The RNG algorithm is pinned along with the seed. BiocParallel switches
  # its workers to the L'Ecuyer-CMRG generator, so a bare set.seed() inside a
  # worker would give a different random stream than the same call in a
  # serial run (and R < 3.6 used a different sample.kind), making segmentation
  # differ between --threads 1 and --threads > 1 or between R versions.
  sample_seed <- opt$seed + sum(utf8ToInt(sample_name))
  set.seed(
    sample_seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  msg("Processing sample: ", sample_name)
  msg("  Seed for this sample: ", sample_seed, " (DNAcopy CBS uses permutation testing)")
  msg("  BAM: ", bamfile)

  include_Y <- sex_map$sex[match(sample_name, sex_map$sample)] == "XY"

  chrom_filter <- if (include_Y) NA else "Y"

  rc <- binReadCounts(
    bins = bins,
    bamfiles = bamfile,
    bamnames = sample_name,
    minMapq = opt$min_mapq,
    verbose = opt$verbose
  )

  # Reads (MAPQ >= --min-mapq) counted into genomic bins, before any bin
  # filtering. Equals QDNAseq's pData(rc)$total.reads (colSums of counts).
  total_reads <- sum(Biobase::assayDataElement(rc, "counts"), na.rm = TRUE)

  rc <- applyFilters(
    rc,
    residual = TRUE,
    blacklist = TRUE,
    mappability = opt$map_cutoff,
    bases = opt$bases_cutoff,
    chromosomes = chrom_filter,
    verbose = opt$verbose
  )

  zero_bins <- Biobase::assayDataElement(rc, "counts")[, 1] == 0
  n_zero <- sum(zero_bins & Biobase::fData(rc)$use)
  Biobase::fData(rc)$use <- Biobase::fData(rc)$use & !zero_bins
  if (n_zero > 0) {
    msg("  Excluded ", n_zero, " zero-read bins that otherwise passed filters")
  }
  
  cn <- estimateCorrection(rc)
  cn <- correctBins(cn)
  cn <- normalizeBins(cn)
  cn <- smoothOutlierBins(cn)
  cn <- segmentBins(
    cn,
    alpha = opt$segment_alpha,
    undo.splits = "sdundo",
    undo.SD = opt$segment_undo_sd,
    verbose = opt$verbose
  )

  apply_ploidy_correction <- include_Y && !isTRUE(opt$no_sex_ploidy_correction)

  plot_df <- extract_plot_df(
    cn,
    include_Y = include_Y,
    is_xy = apply_ploidy_correction
  ) |>
    filter(sample == sample_name) |>
    flag_excluded_bins(exclude_regions) |>
    add_qc_flags(
      bases_cut = opt$bases_cutoff,
      map_cut = opt$map_cutoff
    )

  offs <- get_chr_offsets(cn, include_Y = include_Y)

  plot_df_pos <- add_pos_cum(plot_df, offs)

  ok_bins <- plot_df_pos |>
    filter(
      status == "OK",
      is.finite(CN_log2),
      is.finite(pos_cum)
    )

  filtered_bins <- plot_df_pos |>
    filter(
      status != "OK",
      is.finite(CN_log2),
      is.finite(pos_cum)
    )

  seg_input <- plot_df_pos |>
    filter(
      is.finite(SEG_log2),
      is.finite(pos_cum),
      is.finite(pos_cum_end)
    )

  seg_runs <- compute_segment_runs(seg_input, tol = opt$segment_tol)

  p <- plot_genomewide(
    ok_df = ok_bins,
    filtered_df = filtered_bins,
    seg_runs = seg_runs,
    offs = offs,
    ymin = opt$ymin,
    ymax = opt$ymax,
    plot_filtered_bins = isTRUE(opt$plot_filtered_bins)
  )

  out_plot <- file.path(opt$out_dir, paste0(sample_name, ".plot.png"))
  png_type <- if (isTRUE(capabilities("aqua"))) "quartz" else "cairo"
  grDevices::png(
    filename = out_plot,
    width = opt$width_mm / 25.4,
    height = opt$height_mm / 25.4,
    units = "in",
    res = 600,
    bg = "white",
    type = png_type
  )
  print(p)
  grDevices::dev.off()

  out_rds <- file.path(opt$out_dir, paste0(sample_name, ".qdnaseq.rds"))
  saveRDS(cn, out_rds)

  out_bins <- file.path(opt$out_dir, paste0(sample_name, ".bins.tsv.gz"))
  export_table_gz(
    plot_df_pos |>
      mutate(start = as.integer(start), end = as.integer(end)) |>
      select(
        sample, chromosome, start, end, bases_pct, mappability_pct,
        low_bases, low_map, excluded_region, status, CN_raw, SEG_raw,
        CN_log2, SEG_log2, pos_cum, pos_cum_end
      ),
    out_bins
  )

  out_segments <- file.path(opt$out_dir, paste0(sample_name, ".segments.tsv"))
  seg_out <- seg_runs |>
    transmute(
      sample      = sample_name,
      chromosome,
      start       = as.integer(start),
      end         = as.integer(end),
      size_mb     = round((end - start + 1) / 1e6, 2),
      n_bins,
      log2_ratio  = round(seg_y, 3),
      copy_ratio  = round(2^seg_y, 3)
    )
  write.table(seg_out, file = out_segments, sep = "\t", row.names = FALSE, quote = FALSE)

  qc_row <- sample_qc_summary(
    plot_df_pos,
    total_reads = total_reads,
    seg_runs = seg_runs
  ) |>
    mutate(
      sample = sample_name,
      sex_mode = if (include_Y) "XY" else "XX",
      sex_ploidy_corrected = apply_ploidy_correction,
      plot_file = basename(out_plot),
      qdnaseq_rds = basename(out_rds),
      bins_tsv_gz = basename(out_bins),
      segments_tsv = basename(out_segments)
    )

  msg("  Wrote plot: ", out_plot)
  msg("  Wrote object: ", out_rds)
  msg("  Wrote bins: ", out_bins)
  msg("  Wrote segments: ", out_segments)

  list(qc = qc_row, plot = p)
}

# ---------------------------------------------------------------------------
# Dispatch: serial (default) or parallel across samples
# ---------------------------------------------------------------------------

n_workers <- 1L

if (!is.null(opt$threads) && opt$threads > 1) {
  if (.Platform$OS.type == "windows") {
    msg(
      "--threads ", opt$threads, " was requested, but this script's ",
      "parallel path only supports fork-based workers (Mac/Linux), which ",
      "are not available on Windows. Falling back to serial execution ",
      "(--threads 1) rather than risk an unverified Windows multi-process ",
      "code path producing silently incorrect results."
    )
  } else {
    n_workers <- max(1L, min(
      as.integer(opt$threads),
      length(bamfiles),
      max(1L, parallel::detectCores() - 1L)
    ))
  }
}

if (n_workers > 1) {
  msg(
    "Running ", length(bamfiles), " sample(s) across ", n_workers,
    " parallel worker(s) (MulticoreParam). Per-sample progress messages ",
    "below may interleave between workers; this does not affect correctness."
  )
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers)
  qc_rows <- BiocParallel::bplapply(
    seq_along(bamfiles),
    process_one_sample,
    bamfiles = bamfiles,
    bamnames = bamnames,
    sex_map = sex_map,
    bins = bins,
    exclude_regions = exclude_regions,
    opt = opt,
    BPPARAM = BPPARAM
  )
} else {
  qc_rows <- lapply(
    seq_along(bamfiles),
    process_one_sample,
    bamfiles = bamfiles,
    bamnames = bamnames,
    sex_map = sex_map,
    bins = bins,
    exclude_regions = exclude_regions,
    opt = opt
  )
}

qc_summary <- bind_rows(lapply(qc_rows, `[[`, "qc")) |>
  select(sample, everything())

plots <- stats::setNames(lapply(qc_rows, `[[`, "plot"), bamnames)
summary_pdf <- file.path(opt$out_dir, "CNAscreen_summary.pdf")
write_run_summary_pdf(plots, summary_pdf)
msg("Wrote run summary: ", summary_pdf)

# Write read counts as plain integers (write.table would otherwise print
# round values such as 50000000 as 5e+07).
qc_summary$total_reads <- format(
  qc_summary$total_reads,
  scientific = FALSE,
  trim = TRUE
)

write.table(
  qc_summary,
  file = file.path(opt$out_dir, "sample_qc_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write_session_info(opt$out_dir, opt$seed)

msg("Done.")

