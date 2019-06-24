#' @include Ensemble.SDM.R Stacked.SDM.R checkargs.R
#' @importFrom raster writeRaster
NULL

#' Save ensemble SDMs and SSDMs
#' \strong{(Esp)} Guarda SDMs ensamblados y SSDMs
#'
#' Allows to save S4 \linkS4class{Ensemble.SDM} and \linkS4class{Stacked.SDM}
#' class objects.
#' \strong{(Esp)} Permite guardar objetos de clase S4 \linkS4class{Ensemble.SDM} y
#' \linkS4class{Stacked.SDM}.
#'
#' @param enm Ensemble.SDM. Ensemble SDM to be saved.
#' \strong{(Esp)} SDM ensambladoa ser giardado.
#' @param stack Stacked.SDM. SSDM to be saved.
#' \strong{(Esp)} SSDM a ser guardado.
#' @param name character. Folder name of the model to save.
#' \strong{(Esp)} Nombre de la carpeta donde se guardara el modelo.
#' @param path character. Path to the directory chosen to save the SDM,
#'   by default the path to the current directory.
#' \strong{(Esp)} Ruta del directorio elegido para guardar el SDM, por defecto
#'   la ruta es el actual directorio.
#' @param verbose logical. If set to true, allows the function to print text in the
#'   console.
#' \strong{(Esp)} Si es verdadero (\code{TRUE}) permite la función de escribir texto
#'   en la consola.
#' @param GUI logical. Don't take that argument into account (parameter for the
#'   user interface).
#' \strong{(Esp)} No tomar en cuenta este argumento (parámetro para la interfaz de usuario)
#'
#' @return Nothing in R environment. Creates folders, tables and rasters
#'   associated to the SDM. Tables are in .csv and rasters in .grd/.gri.
#' \strong{(Esp)} Nada en el ambiente (environment) de R. Crea una carpeta, tablas y rásters
#'   asociados al SDM. TAblas estan en formato .csv y raster en .grd/.gri
#'
#' @seealso \code{\link{load.model}}
#'
#' @name save.model
#'
NULL

#' @rdname save.model
#' @export
setGeneric('save.enm', function (enm, name = strsplit(enm@name, '.', fixed = TRUE)[[1]][1],
                                 path = getwd(), verbose = TRUE, GUI = FALSE) {return(standardGeneric('save.enm'))})

#' @rdname save.model
#' @export
setMethod('save.enm', 'Ensemble.SDM', function (enm,
                                                        name = strsplit(enm@name, '.Ensemble.SDM', fixed = TRUE)[[1]][1],
                                                        path = getwd(),
                                                        verbose = TRUE, GUI = FALSE) {
  # Check arguments
  .checkargs(enm = enm, name = name, path = path, verbose = verbose, GUI = GUI)

  if (verbose) {
    cat("Saving ensemble model results \n")
  }
  # Directories creation
  dir.create(path = paste0(path, "/", name))
  dir.create(path = paste0(path, "/", name, "/Rasters"))
  dir.create(path = paste0(path, "/", name, "/Tables"))

  # Raster saving
  if (verbose) {
    cat("   rasters ...")
  }
  writeRaster(enm@projection[[1]], paste0(path, "/", name, "/Rasters/Probability"),
              "GTiff", overwrite = TRUE)
  writeRaster(enm@binary[[1]], paste0(path, "/", name, "/Rasters/Binary"),
              "GTiff", overwrite = TRUE)
  writeRaster(enm@uncertainty, paste0(path, "/", name, "/Rasters/uncertainty"),
              "GTiff", overwrite = TRUE)
  if (verbose) {
    cat("saved \n")
  }

  # Tables saving
  if (verbose) {
    cat("   tables ...")
  }
  write.csv(enm@evaluation, paste0(path, "/", name, "/Tables/ENMeval.csv"))
  write.csv(enm@algorithm.evaluation, paste0(path, "/", name, "/Tables/AlgoEval.csv"))
  write.csv(enm@algorithm.correlation, paste0(path, "/", name, "/Tables/AlgoCorr.csv"))
  write.csv(enm@variable.importance, paste0(path, "/", name, "/Tables/VarImp.csv"))
  write.csv(enm@data, paste0(path, "/", name, "/Tables/Data.csv"))
  write.csv(enm@name, paste0(path, "/", name, "/Tables/Name.csv"))
  write.csv(enm@parameters, paste0(path, "/", name, "/Tables/Parameters.csv"))
  if (verbose) {
    cat("saved \n \n")
  }
})

#' @rdname save.model
#' @export
setGeneric('save.stack', function (stack, name = 'Stack', path = getwd(), verbose = TRUE, GUI = FALSE) {return(standardGeneric('save.stack'))})

#' @rdname save.model
#' @export
setMethod('save.stack', 'Stacked.SDM', function (stack, name = 'Stack',
                                                                        path = getwd(),
                                                                        verbose = TRUE, GUI = FALSE) {
  # Check arguments
  .checkargs(stack = stack, name = name, path = path, verbose = verbose,
             GUI = GUI)

  if (verbose) {
    cat("Saving stack species model results \n")
  }
  # Directories creation
  dir.create(path = paste0(path, "/", name))
  path = paste0(path, "/", name)
  dir.create(path = paste0(path, "/", "Stack"))
  dir.create(path = paste0(path, "/", "Species"))
  dir.create(path = paste0(path, "/", "Stack", "/Rasters"))
  dir.create(path = paste0(path, "/", "Stack", "/Tables"))

  # Raster saving
  if (verbose) {
    cat("   rasters ...")
  }
  writeRaster(stack@diversity.map, paste0(path, "/", "Stack", "/Rasters/Diversity"),
              "GTiff", overwrite = TRUE)
  writeRaster(stack@endemism.map, paste0(path, "/", "Stack", "/Rasters/Endemism"),
              "GTiff", overwrite = TRUE)
  writeRaster(stack@uncertainty, paste0(path, "/", "Stack", "/Rasters/uncertainty"),
              "GTiff", overwrite = TRUE)
  cat("saved \n")

  # Tables saving
  if (verbose) {
    cat("   tables ...")
  }
  write.csv(stack@evaluation, paste0(path, "/", "Stack", "/Tables/StackEval.csv"))
  write.csv(stack@algorithm.evaluation, paste0(path, "/", "Stack", "/Tables/AlgoEval.csv"))
  write.csv(stack@algorithm.correlation, paste0(path, "/", "Stack", "/Tables/AlgoCorr.csv"))
  write.csv(stack@variable.importance, paste0(path, "/", "Stack", "/Tables/VarImp.csv"))
  write.csv(stack@name, paste0(path, "/", "Stack", "/Tables/Name.csv"))
  write.csv(stack@parameters, paste0(path, "/", "Stack", "/Tables/Parameters.csv"))
  cat("saved \n\n")

  # ENMS saving
  if (verbose) {
    cat("   enms ... \n\n")
  }
  for (i in seq_len(length(stack@enms))) {
    save.enm(stack@enms[[i]], path = paste0(path, "/", "Species"), verbose = verbose,
             GUI = GUI)
  }
  if (verbose) {
    cat("saved \n \n")
  }
})
