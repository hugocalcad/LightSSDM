#' @include Algorithm.SDM.R
#' @import methods
#' @import SpaDES.tools
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

setGeneric("project", function(obj, Env, folder_tmp = NULL, ...) {
  return(standardGeneric("project"))
})

setMethod("project", "Algorithm.SDM", function(obj, Env, folder_tmp = NULL, ...) {
  model = get_model(obj, ...)
  path <- get("tmpdir", envir = .PkgEnv)
  if(!is.null(folder_tmp))
    path = paste0(path, "/.raster/", folder_tmp, "/.split")
  else
    path = paste0(path,"/.raster/.split")
  if(!(file.path(path, names(Env[[1]])) %in% list.dirs(path))){
    for(i in seq_len(length(Env@layers))){
      splitRaster(Env[[i]],1,1,path = paste0(path,"/", names(Env[[i]])))
    }
  }
  n = length(list.files(paste0(path,"/",names(Env[[1]]))))/2
  proj_2 = list()
  bin = list()

  for(i in seq_len(n)){
    listado_completo <- list.files(path , pattern = paste("*tile",i,".gri$", sep = ''),
                                   recursive = TRUE, full.names = TRUE )
    Env_temp <- raster::stack(listado_completo)
    proj_2[[i]] = suppressWarnings(raster::predict(Env_temp, model, fun = function(model, x) {
      x = as.data.frame(x)
      for (i in seq_len(length(Env_temp@layers))) {
        if (Env_temp[[i]]@data@isfactor) {
          x[, i] = as.factor(x[, i])
          x[, i] = droplevels(x[, i])
          levels(x[, i]) = Env_temp[[i]]@data@attributes[[1]]$ID
        }
      }
      return(predict(model, x))
    }))
    # Rescaling projection
    proj_2[[i]] = reclassify(proj_2[[i]], c(-Inf, 0, 0))
    if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
      bin[[i]] = reclassify(proj_2[[i]], c(-Inf,obj@evaluation$threshold,0,
                                       obj@evaluation$threshold,Inf,1))
  }
  if(n==1){
    proj <- proj_2[[n]]
    obj@binary <- bin[[n]]
  }else{
    proj <- mergeRaster(proj_2)
    obj@binary <- mergeRaster(bin)
  }
  if(all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj

  return(obj)
})

setMethod("project", "MAXENT.SDM", function(obj, Env, folder_tmp = NULL, ...) {
  model = get_model(obj, Env, ...)
  path <- get("tmpdir", envir = .PkgEnv)
  if(!is.null(folder_tmp))
    path = paste0(path, "/.raster/", folder_tmp, ".split")
  else
    path = paste0(path,"/.raster/.split")
  if(!(file.path(path, names(Env[[1]])) %in% list.dirs(path))){
    for(i in seq_len(length(Env@layers))){
      splitRaster(Env[[i]],1,1,path = paste0(path,"/", names(Env[[i]])))
    }
  }
  n = length(list.files(paste0(path,"/",names(Env[[1]]))))/2
  proj_2 = list()
  bin = list()

  for(i in seq_len(n)){
    listado_completo <- list.files(path , pattern = paste("*tile",i,".gri$", sep = ''),
                                   recursive = TRUE, full.names = TRUE )
    Env_temp <- raster::stack(listado_completo)
    proj_2[[i]] = raster::predict(Env_temp, model, fun = function(model, x) {
      x = as.data.frame(x)
      for (i in seq_len(length(Env_temp@layers))) {
        if (Env_temp[[i]]@data@isfactor) {
          x[, i] = as.factor(x[, i])
          x[, i] = droplevels(x[, i])
          levels(x[, i]) = Env_temp[[i]]@data@attributes[[1]]$ID
        }
      }
      return(predict(model, x))
    })
    # Rescaling projection
    proj_2[[i]] = reclassify(proj_2[[i]], c(-Inf, 0, 0))
    if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
      bin[[i]] <- reclassify(proj_2[[i]], c(-Inf,obj@evaluation$threshold,0,
                                       obj@evaluation$threshold,Inf,1))
  }

  if(n==1){
    proj <- proj_2[[n]]
    obj@binary <- bin[[i]]
  }else{
    proj <- mergeRaster(proj_2)
    obj@binary <- mergeRaster(bin)
  }

  if(!all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj

  return(obj)
})
