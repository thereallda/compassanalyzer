#' @rdname createCompass
#' @export
#'
setClass(
  Class = "Compass",
  contains = "SummarizedExperiment",
  slots = list(
    counts = "list",
    norm_factors = "list",
    ratio = "list",
    ratio_filtered = "list",
    parameter = "list"
  )
)

#' Compass object and constructor
#'
#' @description \code{Compass} object extends the \code{SummarizedExperiment} class.
#' The \code{createCompass} is a easy constructor of \code{Compass} object
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param bio.group Vector of samples group, e.g., c("Young.Input","Young.Enrich","Old.Input","Old.Enrich").
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param batch.group Vector of samples batch, e.g., c("A","A","B","B"), default: NULL.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g., "^FB" stands for gene id of fly spike-in , default: NULL.
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#' @param synthetic.id Vector of synthetic RNA id, e.g. c("Syn1","Syn2"), default: NULL.
#'
#' @return Compass object
#'
#' @details Description of each slot:
#'   \code{assay} \code{SummarizedExperiment::Assays} object, contains all counts.
#'   \code{counts} list for holding raw/normalized counts of sample and spike_in.
#'   \code{norm_factors} list for holding normalization factors of sample and spike_in.
#'   \code{ratio} list for holding the ratio of all genes.
#'   \code{ratio_filtered} lists for holding the ratio of genes with ratios smaller than 1 in all samples.
#'   \code{parameter} list of parameters.
#'
#' @export
#'
#' @importFrom methods validObject
#' @importFrom enONE countReplicate
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
createCompass <- function(data,
                        bio.group,
                        enrich.group,
                        batch.group = NULL,
                        spike.in.prefix = NULL,
                        input.id = "Input",
                        enrich.id = "Enrich",
                        synthetic.id = NULL
                        ) {

  # parameters
  params <- list(
    spike.in.prefix = spike.in.prefix,
    synthetic.id = synthetic.id,
    input.id = input.id,
    enrich.id = enrich.id
  )

  # create Assay object
  ## count_assay is a list for holding raw/normalized counts of sample and spike_in
  count_assay <- list(sample = list(), spike_in = list())
  ## ratio_assay and ratio_filtered_assay are data.frame for holding all/filtered ratio
  ratio_assay <- ratio_filtered_assay <- list(sample=data.frame(), spike_in=data.frame())

  # norm_factors slot
  norm_factors <- list(sample = list(), spike_in = list())

  # SummarizedExperiment object
  ## rowData for mapping gene id
  rowDf <- S4Vectors::DataFrame(GeneID=rownames(data),
                                SpikeIn=(rownames(data) %in% grep(spike.in.prefix,rownames(data),value=TRUE))
  )
  # if synthetic RNA id provided
  if (!is.null(synthetic.id)) {
    rowDf$Synthetic <- rowDf$GeneID %in% grep(paste(synthetic.id,collapse = "|"),rowDf$GeneID,value=TRUE)
  } else {
    rowDf$Synthetic <- rep(FALSE, nrow(rowDf))
  }
  ## colData for mapping samples
  colDf <- S4Vectors::DataFrame(
    id = colnames(data),
    condition = bio.group,
    enrich = enrich.group,
    replicate = enONE::countReplicate(bio.group),
    batch = NA_character_
  )

  if (!is.null(batch.group)) colDf$batch <- batch.group
  # create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(as.matrix(data)),
                                                   colData = colDf,
                                                   rowData = rowDf)
  # create Compass object
  Compass <- methods::new("Compass",
                       se,
                       counts = count_assay,
                       norm_factors = norm_factors,
                       ratio = ratio_assay,
                       ratio_filtered = ratio_filtered_assay,
                       parameter = params
                       )

  # separate the raw counts
  Compass@counts$sample[["Raw"]] <- data[grep(paste(c(spike.in.prefix, synthetic.id), collapse = "|"), rownames(data), invert = TRUE),]
  Compass@counts$spike_in[["Raw"]] <- data[grep(spike.in.prefix, rownames(data)),]

  validObject(Compass)
  return(Compass)
}

setValidity("Compass", function(object) {

  # check input.id and enrich.id match the enrich.group
  if (!all(c(object@parameter$input.id, object@parameter$enrich.id) %in% object$enrich)) {
    return("The `input.id` and/or `enrich.id` must be the same as `enrich.group`.")
  }

  # check spike.in.prefix match the rownames
  if (length(grep(object@parameter$spike.in.prefix, rownames(object))) < 1) {
    return("The `spike.in.prefix` does not match the rownames in the count matrix.")
  }

  # check synthetic.id match the rownames
  if (!is.null(object@parameter$synthetic.id) & !all(object@parameter$synthetic.id %in% rownames(object))) {
    return("The `synthetic.id` are not presented in the rownames of count matrix.")
  }

  # check length of bio.group match the sample size
  if (length(object$condition) != ncol(object)) {
    return("The number of elements in `bio.group` does not match the sample size.")
  }

  # check length of enrich.group match the sample size
  if (length(object$enrich) != ncol(object)) {
    return("The number of elements in enrich.group` does not match the sample size.")
  }

  # check length of batch.group match the sample size
  if (!all(is.na(object$batch)) & length(object$batch) != ncol(object)) {
    return("The number of elements in `batch.group` does not match the sample size.")
  }
})
