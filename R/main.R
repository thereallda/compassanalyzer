#' Normalization and ratio assessment
#'
#' @param object Compass object
#' @param adjust Whether to perform linear regression-based adjustment, default: TRUE
#' @param prop.top.enrich Proportion of top-enriched genes to use for
#' adjustment of non-specific enrichment, by default all genes are used.
#' @param top.down Whether using top-to-down or down-to-top enriched genes for
#' the adjustment of non-specific enrichment, default: TRUE (top-down)
#' @param pseudo.count A numeric scalar of pseudo-counts to be added to each gene, default: 1.
#'
#' @return Compass object
#' @export
#'
#' @importFrom edgeR calcNormFactors
CompassAnalyze <- function(object,
                           adjust = TRUE,
                           prop.top.enrich = 1,
                           top.down = TRUE,
                           pseudo.count = 1
                          ) {
  # retrieve parameters from Compass object
  bio.group <- object$condition
  enrich.group <- object$enrich

  spike.in.prefix <- object@parameter$spike.in.prefix
  input.id <- object@parameter$input.id
  enrich.id <- object@parameter$enrich.id
  synthetic.id <- object@parameter$synthetic.id

  # get counts
  data <- SummarizedExperiment::assay(object)
  # sample counts
  counts_sam <- data[grep(paste(c(spike.in.prefix, synthetic.id), collapse = "|"), rownames(data), invert = TRUE),]
  # spike-in counts
  counts_sp <- data[grep(spike.in.prefix, rownames(data)),]

  # Global scaling
  cat("Global scaling...\n")
  scale_factors <- calcScaleFactor(data,
                                   spike.in.prefix = spike.in.prefix,
                                   enrich.group = enrich.group,
                                   input.id = input.id,
                                   enrich.id = enrich.id
                                   )
  counts_sam_scale <- t(t(counts_sam) * scale_factors)

  # Linear regression-based adjustment
  if (adjust) {
    cat("Linear regression-based adjustment...\n")
    adjust_factors <- calcAdjustFactor(data,
                                       spike.in.prefix = spike.in.prefix,
                                       enrich.group = enrich.group,
                                       input.id = input.id,
                                       enrich.id = enrich.id,
                                       scale.factor = scale_factors,
                                       prop.top.enrich = prop.top.enrich,
                                       top.down = top.down,
                                       pseudo.count = pseudo.count)
    counts_sam_adjust <- t(t(counts_sam_scale)*adjust_factors)

  } else {
    # no adjustment
    adjust_factors <- rep(1, ncol(data))
    counts_sam_adjust <- NULL
  }

  # Computation of NCIN Ratio
  cat("Computation of NCIN Ratio...\n")
  ratio_ls <- calcNCIN(data,
                       spike.in.prefix = spike.in.prefix,
                       enrich.group = enrich.group,
                       input.id = input.id,
                       enrich.id = enrich.id,
                       scale.factor = scale_factors,
                       adjust.factor = adjust_factors,
                       pseudo.count = pseudo.count
  )

  # save results in object
  object@norm_factors[["sample"]] <- list("scaled" = scale_factors,
                                          "adjusted" = adjust_factors)

  object@ratio[["sample"]] <- ratio_ls$ratio
  object@ratio_filtered[["sample"]] <- ratio_ls$ratio_filtered

  # save run parameter in object
  parameter.run <- list(
    adjust = adjust,
    prop.top.enrich = prop.top.enrich,
    top.down = top.down,
    pseudo.count = pseudo.count
  )
  object@parameter <- c(object@parameter, parameter.run)

  # save counts in object
  # store 'Raw' count matrix
  Counts(object, slot = "sample", method = "Raw") <- counts_sam
  Counts(object, slot = "spike_in", method = "Raw") <- counts_sp

  # store normalized count matrix
  Counts(object, slot = "sample", method = "scaled") <- counts_sam_scale
  Counts(object, slot = "sample", method = "adjusted") <- counts_sam_adjust

  validObject(object)
  return(object)
}

#' Generate scale factors based on spike-ins
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g.,
#'  "^FB" stands for gene id of fly spike-in, default: NULL.
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#'
#' @return Vector of scaling factors
#' @export
#'
#' @importFrom edgeR calcNormFactors
calcScaleFactor <- function(data,
                            spike.in.prefix = NULL,
                            enrich.group,
                            input.id = "Input",
                            enrich.id = "Enrich") {

  # initialize parameters
  input.idx <- grep(input.id, enrich.group)
  enrich.idx <- grep(enrich.id, enrich.group)

  # get counts from spike-in
  counts_spk <- data[grep(spike.in.prefix, rownames(data)),]

  # get counts from sample
  counts_sam <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]

  # compute scale factors using TMM method
  nf1 <- edgeR::calcNormFactors(counts_spk, refColumn = input.idx)
  sf1 <- nf1/colSums(counts_spk)*1e6

  # adjust scale factor by ratio of spike-in/non-spike-in
  lib_ratio <- colSums(counts_spk)/colSums(data)
  lib_ratio <- lib_ratio/(1-lib_ratio)
  lib_ratio <- c(lib_ratio[input.idx]/mean(lib_ratio[input.idx]),
                 lib_ratio[enrich.idx]/mean(lib_ratio[enrich.idx]))
  scale.factor <- sf1*lib_ratio

  return(scale.factor)
}

#' Generate adjust factors based on data after scaling
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g.,
#'  "^FB" stands for gene id of fly spike-in, default: NULL.
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#' @param scale.factor Vector of scale factors generated from \code{calcScaleFactor}.
#' @param prop.top.enrich Proportion of top-enriched genes to use for
#' adjustment of non-specific enrichment, by default all genes are used.
#' @param top.down Whether using top-to-down or down-to-top enriched genes for
#' the adjustment of non-specific enrichment, default: TRUE (top-down)
#' @param pseudo.count A numeric scalar of pseudo-counts to be added to each gene, default: 1.
#'
#' @return Vector of adjust factors
#' @export
#'
#' @importFrom stats lm
#' @importFrom utils head
calcAdjustFactor <- function(data,
                             spike.in.prefix = NULL,
                             enrich.group,
                             input.id = "Input",
                             enrich.id = "Enrich",
                             scale.factor,
                             prop.top.enrich = 1,
                             top.down = TRUE,
                             pseudo.count = 1) {
  # initialize parameters
  input.idx <- grep(input.id, enrich.group)
  enrich.idx <- grep(enrich.id, enrich.group)

  # get counts from sample
  counts_sam <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]

  # # ensure each samples have at least one count
  # counts_sam <- counts_sam[rowSums(counts_sam > 0) == ncol(counts_sam),]

  # add pseudo counts to shrink ratio
  counts_sam <- addPseudoCount(counts_sam, pseudo.count = pseudo.count)

  # scale sample counts by scale factor
  counts_sam_scale <- t(t(counts_sam)*scale.factor)

  # calculate empirical NCIN ratio
  eratio_df <- counts_sam_scale[,enrich.idx]/counts_sam_scale[,input.idx]

  # get all non-specific enriched genes based on empirical ratio
  neg_all <- eratio_df[rowSums(eratio_df>1) > 0.7*ncol(eratio_df),]
  neg_avg <- rowMeans(neg_all)
  neg_all <- cbind(neg_all, neg_avg)

  # order neg based on their average in decreasing
  neg_all <- neg_all[order(neg_all[,"neg_avg"], decreasing = top.down),]

  # used prop.top.enrich of non-specific enriched genes for adjustment
  neg_all <- utils::head(neg_all, n = floor(nrow(neg_all)*prop.top.enrich))

  # create centroid based on non-specific enriched genes from input samples
  centroid <- rowMeans(counts_sam_scale[rownames(neg_all), input.idx])

  # linear regression between centroid and each sample
  reg_df <- log2(as.data.frame(cbind(counts_sam_scale[rownames(neg_all),], centroid)))
  adjust.factor <- c()
  for (id in colnames(counts_sam)) {
    dfi <- reg_df[,c("centroid", id)]
    colnames(dfi) <- c("centroid", "x")
    lmi <- stats::lm(centroid ~ 1 + offset(x), data = dfi)
    adjust.factor[id] <- 2^lmi$coefficients
  }

  return(adjust.factor)
}

#' Compute NCIN ratio
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g.,
#'  "^FB" stands for fly spike-in id, default: NULL.
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#' @param scale.factor Vector of scale factors generated from \code{calcScaleFactor}.
#' @param adjust.factor Vector of adjust factors generated from \code{calcAdjustFactor}.
#' @param pseudo.count A numeric scalar of pseudo-counts to be added to each gene, default: 1.
#'
#' @return List containing the unfiltered and filtered ratio of each gene from each sample.
#' @export
#'
calcNCIN <- function(data,
                     spike.in.prefix = NULL,
                     enrich.group,
                     input.id = "Input",
                     enrich.id = "Enrich",
                     scale.factor,
                     adjust.factor,
                     pseudo.count = 1) {

  # initialize parameters
  input.idx <- grep(input.id, enrich.group)
  enrich.idx <- grep(enrich.id, enrich.group)

  # get counts from sample
  counts_sam <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]

  # # ensure each samples have at least one count
  # counts_sam <- counts_sam[rowSums(counts_sam > 0) == ncol(counts_sam),]

  # add pseudo counts to shrink ratio
  counts_sam <- addPseudoCount(counts_sam, pseudo.count = pseudo.count)

  # scale sample counts by scale factor
  counts_sam_scale <- t(t(counts_sam)*scale.factor)

  # further adjust sample counts by adjust factor
  counts_sam_adj <- t(t(counts_sam_scale)*adjust.factor)

  # calculate quantitative ratio
  ratio_df <- counts_sam_adj[,enrich.idx]/counts_sam_adj[,input.idx]

  # keep gene with ratio smaller than 1 in all samples
  ratio_filtered_df <- ratio_df[rowSums(ratio_df<1) == ncol(ratio_df), ]

  # ratio list
  ratio_ls <- list("ratio" = ratio_df, "ratio_filtered" = ratio_filtered_df)

  return(ratio_ls)
}
