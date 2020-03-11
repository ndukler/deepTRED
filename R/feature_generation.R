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
    (GenomicRanges::end(sites) + GenomicRanges::start(sites))/ 2))
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
    site = rep(factor(seq_along(focal_sites)), each = stepsize),
    bin_idx = rep(bin_idx, length(focal_sites))
  )
  gr <- gr[GenomicRanges::start(gr) > 0]
  # Split into GRangesList
  return(GenomicRanges::split(gr, gr$site))
}

collect_features <- function(sites, bigwig_plus, bigwig_minus,
                             bins = 21L, bin_width = 51L) {
# Iterate over seqlevels and collect data one seqlevel at a time, then reshape
# into data.tables where rows are sites and columns are features

}

