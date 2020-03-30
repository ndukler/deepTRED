#' @title Generate GRangesList of windows
#'
#' @description  Generates a set of windows for each element in a
#' \code{\link[GenomicRanges]{GRanges-class}} object. Uses the center of each
#' region as the 0th co-ordinate in generating windows.
#'

#' @param sites a \link[GenomicRanges]{GRanges-class} containing the sites of
#' interest. The central site of this range will be the focal site used to
#' generate bins
#' @param bins a single value or vector of length equal to bin_width indicating
#' the number of bins to generate, with the central bin centered on the
#' focal site, must be odd
#' @param bin_width a vector of one or more odd valued integers indicating the
#' width of the bins to be generated.
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name window_regions
#'
#' @export
window_regions <- function(sites, bins = 21L, bin_width = 51L) {
  ##** Begin checks **##
  if (any(bins %% 2 == 0)) {
    stop("bins must be odd-valued")
  }

  if (!class(sites) == "GRanges") {
    stop("sites must be a GRanges object")
  }

  if (!is.integer(bin_width) && bin_width != as.integer(bin_width)) {
    stop("bin_width must be an integer")
  }

  if (any(bin_width %% 2 == 0)) {
    stop("all bin_width values must odd")
  }

  if (length(bins) > 1 && length(bins) != length(bin_width)) {
    stop("If multiple values are specified for bins there must be an equal ",
         "number of bin_width values specified")
  }

  ##** End checks **##

  # Make sure length of bins matches length of bin_width
  if (length(bins) != length(bin_width)) {
    bins <- rep(bins, length(bin_width))
  }
  # Compute focal sites
  focal_sites <- as.integer(round(
    (GenomicRanges::end(sites) + GenomicRanges::start(sites)) / 2))
  # Preallocate start vector
  bin_starts <- integer(sum(bins * length(sites)))
  # Pre-compute start adjustments and unlist
  start_adjust <- unlist(mapply(function(bins, bin_width) {
    as.integer(seq(from = -(bin_width - 1) / 2 - ((bins - 1) / 2) * bin_width,
        by = bin_width, length.out = bins))
  }, bins = as.list(bins), bin_width = as.list(bin_width), SIMPLIFY = FALSE))
  # Generate vector of widths per site
  widths <- rep(bin_width, times = bins)
  # Pre-calculate the bin column indicies
  bin_idx = unlist(lapply(as.list(bins), function(x) seq_len(x)))
  # Iterate over start sites
  idx <- 1
  stepsize <- sum(bins)
  for(i in seq_along(focal_sites)) {
    bin_starts[idx:(idx + stepsize - 1)] <-
      as.integer(focal_sites[i] + start_adjust)
    idx <- idx + stepsize
  }
  # Create GRanges object out of the starts
  gr <- GenomicRanges::GRanges(
    seqnames = rep(GenomicRanges::seqnames(sites), each = stepsize),
    ranges = IRanges::IRanges(start = bin_starts,
                              width = rep(widths, length(focal_sites))),
    strand = rep(GenomicRanges::strand(sites), each = stepsize),
    site = rep(factor(seq_along(focal_sites)), each = stepsize),
    bin_idx = rep(bin_idx, length(focal_sites))
  )
  gr <- gr[GenomicRanges::start(gr) > 0]
  # Split into GRangesList
  return(GenomicRanges::split(gr, gr$site))
}

#' @title Generate GRangesList of windows
#'
#' @description  Generates a set of windows for each transcript in a
#' \code{\link[GenomicRanges]{GRanges-class}} object. Uses the TSS of each
#' region as the 0th co-ordinate in generating windows.
#'
#' @param trancripts a code{\link[GenomicRanges]{GRanges-class}} object
#' containing
#' @param bigwig_plus the path to a bigwig for reads on the plus strand
#' @param bigwig_minus the path to a bigwig for reads on the minus strand
#' @param remove_incomplete remove regions where any window features go outside
#' the chromosomal range
#' @inheritParams window_regions
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name collect_tss_features
#'
#' @export
collect_tss_features <- function(transcripts, bigwig_plus, bigwig_minus,
                             bins = 21L, bin_width = 51L,
                             remove_incomplete = TRUE) {
  # Get TSS of transcripts
  tss <- GenomicRanges::promoters(transcripts, upstream = 0, downstream = 1)
  # Window sites
  windows <- window_regions(tss)
  names(windows) <- as.character(1:length(windows))
  # Filter out incomplete feature vectors
  if (remove_incomplete) {
    repl_bins <- rep(bins, length.out = length(bin_width))
    expected_length <- sum(bins)
    remove_sites <- which(S4Vectors::elementNROWS(windows) != expected_length)
    if(length(remove_sites) > 0) {
      windows <- windows[-remove_sites]
      # Filter down TSS to only the sites that were retained
      tss <- tss[-remove_sites]
    }
    message("Removed ",length(remove_sites),
            " sites with incomplete feature vectors")
  } else {
    stop("Not able to handle ragged input vectors yet")
  }

  if(length(windows) == 0) {
    stop("No sites remaining")
  }

  # Get strand for each feature set
  transcript_strand <- as.vector(GenomicRanges::strand(tss))

  # Collect sense and antisense counts
  # Reverse bin order such that if windows were indexed ..., -1, 0, 1, ... the
  # sense strand always runs in the positive direction and the anti-sense strand
  # always runs in the negative direction
  # S4Vectors::revElements()
  is_sense_positive <- (transcript_strand == '+')
  sense_counts <- c(
    summarize_bigwig(bigwig_plus, windows[is_sense_positive], "sum"),
    summarize_bigwig(bigwig_minus,
                     S4Vectors::revElements(windows[!is_sense_positive]),
                     "sum"))
  sense_counts <-
    do.call("rbind",
            sense_counts[order(as.integer(names(sense_counts)))]
    )
  antisense_counts <- c(
    summarize_bigwig(bigwig_plus,
                     S4Vectors::revElements(windows[!is_sense_positive]),
                     "sum"),
    summarize_bigwig(bigwig_minus, windows[is_sense_positive], "sum"))
  antisense_counts <-
    do.call("rbind",
            antisense_counts[order(as.integer(names(antisense_counts)))]
    )

  # Get row means
  mean_scale <- (matrixStats::rowMeans2(sense_counts) +
                   matrixStats::rowMeans2(antisense_counts)) / 2
  # Remove rows with mean zero as that implies all entries are 0
  keep <- which(mean_scale > 0)
  if (length(keep) != length(mean_scale)) {
    message("Removing sites with all zero entries in their feature vector: ",
            length(keep), " of ", length(mean_scale), " remaining")
    sense_counts <- sense_counts[keep, ]
    antisense_counts <- antisense_counts[keep, ]
    mean_scale <- mean_scale[keep]
    tss <- tss[keep]
  } else {
    stop("No non-zero entries in feature matrix")
  }

  # Get row standard deviations
  sd_scale <- matrixStats::rowSds(cbind(sense_counts, antisense_counts),
                                  center = mean_scale)

  # Scale feature matricies
  plus_counts <- (sense_counts - mean_scale) / sd_scale
  minus_counts <- (antisense_counts - mean_scale) / sd_scale
  attr(sense_counts, "scaled:center") <- mean_scale
  attr(antisense_counts, "scaled:center") <- mean_scale
  attr(sense_counts, "scaled:scale") <- sd_scale
  attr(antisense_counts, "scaled:scale") <- sd_scale
  return(list(sense = sense_counts, antisense = antisense_counts, tss = tss))
}
