#' @rdname Counts
#' @export
setGeneric("Counts", function(object, slot=c("sample","spike_in"), method) standardGeneric("Counts"))

#' @rdname Counts
#' @export
setGeneric("Counts<-", function(object, slot=c("sample","spike_in"), method, value) standardGeneric("Counts<-"))

#' @rdname getFactor
#' @export
setGeneric("getFactor", function(object, slot=c("sample","spike_in"), method) standardGeneric("getFactor"))

#' @rdname getParameter
#' @export
setGeneric("listParameter", function(object) standardGeneric("listParameter"))

#' @rdname getParameter
#' @export
setGeneric("getParameter", function(object, name) standardGeneric("getParameter"))

#' @rdname getRatio
#' @export
setGeneric("getRatio", function(object, slot=c("sample","spike_in"), filter=FALSE) standardGeneric("getRatio"))
