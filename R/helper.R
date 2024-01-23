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

  ratio_long <- ratio.df[syn.meta$id,] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    tidyr::pivot_longer(cols = -id,
                 names_to = "sample",
                 values_to = "ratio") %>%
    dplyr::left_join(syn.meta, by = "id")

  cor_value <- stats::cor.test(ratio_long$ratio, ratio_long$per)
  cor_label <- paste0("R = ", round(cor_value$estimate,3), "; P = ", signif(cor_value$p.value,3))

  # get limits of x and y
  limits.xy <- range(c(ratio_long$per,ratio_long$ratio))

  # expand the range to draw error bar
  limits.xy <- c(limits.xy[1]-0.01, limits.xy[2]+0.01)

  p1 <- ggplot(ratio_long, aes(per, ratio)) +
    stat_summary(geom = "point", fun = mean) +
    stat_summary(fun.data = mean_se, width = 0.008, geom = "errorbar") +
    geom_abline(slope=1) +
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

# For adjusting no visible binding
## synScatter
utils::globalVariables(c("per", "ratio"))
