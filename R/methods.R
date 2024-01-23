#' Accessor for the "counts" slot of Compass object. 
#'
#' @param object Compass. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which counts matrix to get, must be one of the raw or normalized counts matrix presented in the selected slot. 
#' @param value Raw or normalized counts matrix. 
#' @name Counts
#' @aliases Counts Counts,Compass,character,character-method 
#' Counts<-,Compass,character,character,matrix-method
#' 
#' @return matrix. 
#' @export
#'
setMethod("Counts", signature = signature(object="Compass", slot="character", method="character"), 
          function(object, slot=c("sample","spike_in"), method) {
            
            slot <- match.arg(slot, choices = c("sample","spike_in"))
            
            if (is.null(names(object@counts[[slot]]))) {
              stop("Normalizations for ", slot, " not found. At least one normalization should be performed.")
            }
            method <- match.arg(method, choices = names(object@counts[[slot]]))
            
            object@counts[[slot]][[method]]
            })

#' @rdname Counts
#' @name Counts
#' @export "Counts<-"
setReplaceMethod("Counts", signature = signature(object="Compass", slot="character", method="character", value="matrix"),
                 function(object, slot=c("sample","spike_in"), method, value) {
                   slot <- match.arg(slot, choices = c("sample","spike_in"))
                   object@counts[[slot]][[method]] <- value
                   methods::validObject(object)
                   return(object)
                 })


#' Accessor of Compass normalization factors
#'
#' @param object Compass. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods to get, must be one of the methods presented in the selected slot. 
#' @name getFactor
#' @aliases getFactor getFactor,Compass,character,character-method
#'
#' @return vector or list of factors. 
#' @export
#'
setMethod("getFactor", signature = signature(object="Compass", slot="character", method="character"),
          function(object, slot=c("sample","spike_in"), method) {
            slot <- match.arg(slot, choices = c("sample","spike_in"))
            object@norm_factors[[slot]][[method]]
          })


#' Accessor of Compass parameter
#'
#' List all parameters. 
#'
#' @param object Compass. 
#' @name getParameter
#' @aliases listParameter listParameter,Compass-method
#' 
#' @return Vector of parameter names
#' @export
setMethod("listParameter", signature = signature(object="Compass"),
          function(object) {
            names(object@parameter)
          })


#' Accessor of Compass parameter
#'
#' Get specific parameter.
#'  
#' @param object Compass. 
#' @param name Name of the parameter. 
#' @name getParameter
#' @aliases getParameter getParameter,Compass,character-method
#' 
#' @return Vector of parameter 
#' @export
#'
setMethod("getParameter", signature = signature(object="Compass", name="character"),
          function(object, name) {
            object@parameter[[name]]
          })

#' Accessor of ratio
#'
#' Get all or filtered ratio
#'  
#' @param object Compass. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param filter Whether to get the filtered ratio, default FALSE. 
#' @name getRatio
#' @aliases getRatio getRatio,Compass,character-method
#' 
#' @return list of ratio table 
#' @export
#'
setMethod("getRatio", signature = signature(object="Compass", slot="character"),
          function(object, slot=c("sample","spike_in"), filter=FALSE) {
            if (filter) {
              object@ratio_filtered[[slot]]
            } else {
              object@ratio[[slot]]
            }
          })
