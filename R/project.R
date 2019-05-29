#' @include Algorithm.SDM.R
#' @import methods
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

setMethod("project", "Algorithm.SDM", function(obj, Env, ...) {
  model = get_model(obj, ...)

  proj = suppressWarnings(raster::predict(Env, model, fun = function(model, x) {
      x = as.data.frame(x)
      for (j in seq_len(length(Env@layers))) {
        if (Env[[j]]@data@isfactor) {
          x[, j] = as.factor(x[, j])
          x[, j] = droplevels(x[, j])
          levels(x[, j]) = Env[[j]]@data@attributes[[1]]$ID
        }
      }
      if(obj@name == 'MAXNET.SDM')
        return(predict(model, x, type = "logistic"))
      else
        return(predict(model, x))
  }))
  # Rescaling projectio
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
