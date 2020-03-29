#' @title Build labeled data set
#'
#' @description Builds and labels feature vectors for both active and inactive
#' TSS.
#'
#' @param bins A \code{\link[GenomicRanges]{GRanges-class}} object of annotated
#' transcripts
#' @param active_tss a code{\link[GenomicRanges]{GRanges-class}} object
#' of regions annotated as containing a TSS
#' @param bigwig_plus a character vector of paths to bigWig files containing the
#' positions of polymerase on the plus strand
#' @param bigwig_minus a character vector of paths to bigWig files containing
#' the positions of polymerase on the minus strand
#' @param active_tss_distance_cutoff the maximum distance at which an active tss
#' regions is considered to be associated with a gene
#' @param active_tss_distance_cutoff the minimal distance which an annotated tss
#' must be from an active TSS region to be considered inactive
#' @inheritParams window_regions
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#' @include feature_generation.R
#' @name labeled_feature_set
#'
#' @export
labeled_feature_set <- function(labeled_transcripts, bigwig_plus,
                               bigwig_minus, bins = 51L, bin_width = 51L) {

  # Drop inactive seqlevels to avoid spurious warnings
  GenomeInfoDb::seqlevels(labeled_transcripts) <-
    GenomeInfoDb::seqlevelsInUse(labeled_transcripts)

  # Check that active_tss column is present
  if(!"active_tss" %in% colnames(GenomicRanges::mcols(labeled_transcripts))) {
    stop("labeled_transcripts must contain column active_tss")
  }

  if(!is.logical(GenomicRanges::mcols(labeled_transcripts)$active_tss)) {
    stop("active_tss must be column of type logical")
  }

  # Subset transcript to only those not NA
  pre_filter <- length(labeled_transcripts)
  labeled_transcripts <-
    labeled_transcripts[!is.na(labeled_transcripts$active_tss)]
  post_filter <- length(labeled_transcripts)
  message("Removed NA transcripts: ", post_filter, " of ", pre_filter,
          " remaining")

  if(post_filter == 0) {
    stop("No transcripts remaining")
  }

  # Get TSS of transcripts
  tss <- GenomicRanges::promoters(labeled_transcripts,
                                  upstream = 0, downstream = 1)

  # Get features for the postive and negative sets of sites
  active_features <-
    collect_features(tss[tss$active_tss], bigwig_plus, bigwig_minus)
  inactive_features <-
    collect_features(tss[!tss$active_tss], bigwig_plus, bigwig_minus)
  return(list(active_features, inactive_features))
}

#' @title Call active TSS
#'
#' @description Calls active TSS from *-cap based sequencing data
#' @param transcripts foo
#' @param cap_bigwig_plus foo
#' @param cap_bigwig_minus foo
#' @param ctrl_bigwig_plus foo
#' @param ctrl_bigwig_minus foo
#' @param ctrl_cdf_alpha foo
#' @param radius foo
#' @return Return transcripts with an added active_tss column indicating if the
#' tss is active
#'
#' @include bigwig_utils.R
#' @name label_transcripts
#'
#' @export
label_transcripts <- function(transcripts,
                                cap_bigwig_plus, cap_bigwig_minus,
                                ctrl_bigwig_plus, ctrl_bigwig_minus,
                                ctrl_cdf_alpha = 0.01, radius = 150) {

  # Drop inactive seqlevels to avoid spurious warnings
  GenomeInfoDb::seqlevels(transcripts) <-
    GenomeInfoDb::seqlevelsInUse(transcripts)

  # Get TSS of transcripts
  tss <- GenomicRanges::promoters(transcripts, upstream = radius,
                                  downstream = radius)
  # Default tss to being inactive
  tss$active_tss <- NA

  # Compute total reads in each bigwig file pair
  cap_reads <- total_reads(cap_bigwig_plus) + total_reads(cap_bigwig_minus)
  ctrl_reads <- total_reads(ctrl_bigwig_plus) + total_reads(ctrl_bigwig_minus)
  # Convert to units of millions
  cap_rpm <- cap_reads / 1e6
  ctrl_rpm <- ctrl_reads / 1e6

  # Get *-cap counts around promoters
  tss_plus <- GenomicRanges::strand(tss) == "+"
  cap_p <-
    abs(summarize_bigwig(cap_bigwig_plus, bins = tss[tss_plus])) / cap_rpm
  cap_m <-
    abs(summarize_bigwig(cap_bigwig_minus, bins = tss[!tss_plus])) / cap_rpm
  ctrl_p <-
    abs(summarize_bigwig(ctrl_bigwig_plus, bins = tss[tss_plus])) / ctrl_rpm
  ctrl_m <-
    abs(summarize_bigwig(ctrl_bigwig_minus, bins = tss[!tss_plus])) / ctrl_rpm

  # Create an empirical cdf
  ctrl_ecdf <- stats::ecdf(c(ctrl_p, ctrl_m))

  # Label the positive and negative training examples based on the tail
  # probability of the values observed within radius of the tss under the
  # control empirical CDF
  tss[tss_plus][1 - ctrl_ecdf(cap_p) <= ctrl_cdf_alpha]$active_tss <- TRUE
  tss[!tss_plus][1 - ctrl_ecdf(cap_m) <= ctrl_cdf_alpha]$active_tss <- TRUE
  tss[tss_plus][cap_p == 0]$active_tss <- FALSE
  tss[!tss_plus][cap_m == 0]$active_tss <- FALSE

  # Update transcripts with tss info
  transcripts$active_tss <- tss$active_tss

  # Return transcripts GRanges with active_tss column indicating if the tss is
  # active
  return(transcripts)
}
