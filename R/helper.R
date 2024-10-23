#' Scatter plot of synthetic RNA ratios
#'
#' @param ratio.df data.frame containing the NCIN ratios of synthetic spike-ins.
#' @param syn.meta data.frame of metadata of synthetic spike-ins.
#'
#' @return ggplot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom stats cor.test
#' @importFrom scales percent
synScatter <- function(ratio.df, syn.meta) {

  cor_value <- stats::cor.test(ratio.df$ratio, ratio.df$per)
  cor_label <- paste0("R = ", round(cor_value$estimate,3), "; P = ", signif(cor_value$p.value,3))

  # get limits of x and y
  limits.xy <- range(c(ratio.df$per,ratio.df$ratio+ratio.df$ratio.se))

  # expand the range to draw error bar
  limits.xy <- c(limits.xy[1]-0.01, limits.xy[2]+0.01)

  p1 <- ggplot(ratio.df, aes(per, ratio)) +
    geom_point() +
    # se error bar
    geom_errorbar(aes(ymin=ratio-ratio.se,
                      ymax=ratio+ratio.se), width = 0.008) +
    geom_abline(slope=1, linetype='dashed') +
    theme_classic() +
    theme(axis.text = element_text(color="black")) +
    scale_x_continuous(limits = limits.xy, labels = scales::percent) +
    scale_y_continuous(limits = limits.xy, labels = scales::percent) +
    labs(x = "Expected Ratio", y = "Observed Ratio",
         subtitle = cor_label)

  return(p1)
}


#' Add library-size adjusted pseudo-counts
#' @description Add library-size adjusted pseudo-counts adapted from edgeR
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param pseudo.count A numeric scalar of pseudo-counts to be added to each gene, default: 1.
#'
#' @return counts matrix with adding pseudo counts
#' @export
#'
addPseudoCount <- function(data, pseudo.count = 1) {

  # get library size
  lib.size <- colSums(data)

  # library-size adjusted pseudo-counts
  adj.pseudo.count <- pseudo.count*lib.size/mean(lib.size)
  adj.pseudo.count <- matrix(adj.pseudo.count, nrow = nrow(data),
                             ncol = ncol(data), byrow = T) # expand to a matrix

  # add library-size adjusted pseudo-counts
  data.adj <- data + adj.pseudo.count

  return(data.adj)
}


#' Calculation of synthetic spike-in ratio
#'
#' @param object Compass object.
#' @param ratio.shrinkage Whether to perform shrinkage of ratio, default TRUE.
#' @param syn.meta data.frame of metadata of synthetic spike-ins.
#'
#' @return data.frame of synthetic spike-in ratio
#' @export
#'
#' @importFrom dplyr left_join
#' @importFrom stats sd
synRatio <- function(object, syn.meta, ratio.shrinkage = TRUE) {
  if (is.null(object@ratio) | is.null(object@ratio_shrunk)) {
    stop("Pleass run CompassAnalyzer before calculating synthetic spike-in ratios.")
  }
  syn_id <- object@parameter$synthetic.id
  # get prior ratio
  ratio_prior <- object@ratio$sample
  n_samples <- ncol(ratio_prior)
  if (ratio.shrinkage) {
    ratio_shrunk_ls <- ratioShrinkage(ratio_prior, bio.group = NULL)
    ratio_shrunk <- ratio_shrunk_ls[[1]]
    ratio_df <- ratio_shrunk[,c("GeneID","ratio.shrunk","ratio.shrunk.sd")]
    ratio_df$ratio.shrunk.se <- ratio_df$ratio.shrunk.sd/sqrt(n_samples)
    syn_ratio_df <- subset(ratio_df, GeneID %in% syn_id)
    colnames(syn_ratio_df) <- c("GeneID","ratio","ratio.sd","ratio.se")
  } else {
    ratio_sd <- apply(ratio_prior, 1, stats::sd)
    ratio_se <- ratio_sd/sqrt(n_samples)

    ratio_df <- data.frame(
      GeneID = rownames(ratio_prior),
      ratio = rowMeans(ratio_prior),
      ratio.sd = ratio_sd,
      ratio.se = ratio_se
    )
    syn_ratio_df <- subset(ratio_df, GeneID %in% syn_id)
  }

  # merge with syn.meta
  syn_ratio_df <- left_join(syn_ratio_df, syn.meta, by=c('GeneID'='id'))

  return(syn_ratio_df)
}

# For adjusting no visible binding
## synScatter
utils::globalVariables(c("per", "ratio","ratio.se"))
## synRatio
utils::globalVariables(c("GeneID"))
