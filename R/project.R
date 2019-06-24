#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

setMethod("project", "Algorithm.SDM", function(obj, Env, ...) {
  model = get_model(obj, ...)
  ##copying terms on SSDM 0.2.4
  factors <- sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) Env[[i]]@data@attributes[[1]]$ID)
  factors[sapply(factors, is.null)] <- NULL
  names(factors) <- unlist(sapply(seq_len(length(Env@layers)), function(i)
    if(Env[[i]]@data@isfactor) names(Env[[i]])))
  if(length(factors)==0) factors <- NULL
  if (obj@name == "MAXNET.SDM")
    proj = suppressWarnings(raster::predict(Env, model, factors = factors, type = "logistic"))
  else
   proj = suppressWarnings(raster::predict(Env, model,factors = factors) )
  proj= reclassify(proj, c(-Inf, 0, 0))
  if(all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
    obj@binary = reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                    obj@evaluation$threshold,Inf,1))
  return(obj)
})
