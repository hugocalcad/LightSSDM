#' @include checkargs.R
#' @importFrom shiny incProgress
#' @importFrom raster raster stack res extent crop reclassify as.factor resample readAll
NULL

#'Load environmental variables
#'\strong{(Esp)} Cargar variables ambientales
#'
#'Function to load environmental variables in the form of rasters to perform
#'\code{\link{modelling}}, \code{\link{ensemble_modelling}} or \code{\link{stack_modelling}}.
#'\strong{(Esp)} Función para cargar varaibles ambientales en el formato ráster para realizar
#'\code{\link{modelling}}, \code{\link{ensemble_modelling}} ó \code{\link{stack_modelling}}.
#'
#'@param path character. Path to the directory that contains the environmental variables files.
#'  \strong{(Esp)} Ruta del directorio que contiene las variables ambientales
#'@param files character. Files containing the environmental variables If NULL
#'  (default) all files present in the path in the selected format will be loaded.
#'  \strong{(Esp)} Archivos que contienen las variables ambientales si es NULL (por defecto) todos
#'  los archivos en la ruta seleccionada serán cargadas.
#'@param format character. Format of environmental variables files
#'  (including .grd, .tif, .asc, .sdat, .rst, .nc, .tif, .envi, .bil, .img).
#'  \strong{(Esp)} Formato de los archivos de las variables ambientales
#'  (incluye .grd, .tif, .asc, .sdat, .rst, .nc, .tif, .envi, .bil, .img).
#'@param categorical character. Specify whether an environmental variable is a categorical variable.
#'  \strong{(Esp)} Especifica si una variable ambiental es una variable categorica.
#'@param Norm logical. If set to true, normalizes environmental variables between 0 and 1.
#'  \strong{(Esp)}. Si es verdadero (\code{TRUE}), normaliza las variables ambientales entre 0 y 1.
#'@param tmp logical. If set to true, rasters are
#'  read in temporary file avoiding to overload the random access memory. But
#'  beware: if you close R, temporary files will be deleted.
#'  \strong{(Esp)} Si es verdadero (\code{TRUE}), los rásters son almacenados en un archivo temporal
#'  para evitar la sobrecarga de la memoria RAM. Pero ha que tener cuidado: Si se cierra R,
#'  los archivos temporales serán eliminados.
#'@param verbose logical. If set to true, allows the function to print text in the console.
#'  \strong{(Esp)} Si es verdadero (\code{TRUE}) imprime texto en la consola.
#'@param GUI logical. Do not take that argument into account (parameter for the
#'  user interface).
#'  \strong{(Esp)} No tomar en cuenta este argumento (parámetro para la interfaz de usuario)
#'@param folder_tmp character. File name where  environment variables are stored,
#' only if tmp is true.
#' \strong{(Esp)} Nombre de archivo donde  las variables ambientales son almacenadas,
#' solamente si 'tmp' es verdadero (\code{TRUE}).
#'
#'@return A stack containing the environmental rasters (normalized or not).
#'  \strong{(Esp)} Una pila que contiene los ráster ambientales (normalizado o no).
#'
#' @examples
#' \dontrun{
#' load_var(system.file('extdata',  package = 'LightSSDM'))
#' }
#'
#'@seealso \code{\link{load_occ}} to load occurrences.
#' \strong{(Esp)} para cargar ocurrencias.
#'
#'@export
load_var <- function(path = getwd(), files = NULL, format = c(".grd", ".tif",
                                                              ".asc", ".sdat", ".rst", ".nc", ".envi", ".bil", ".img"), categorical = NULL,
                     Norm = TRUE, tmp = TRUE, verbose = TRUE, GUI = FALSE, folder_tmp = NULL) {
  # Check arguments
  .checkargs(path = path, files = files, format = format, categorical = categorical,
             Norm = Norm, tmp = tmp, verbose = verbose, GUI = GUI, folder_tmp = folder_tmp)

  # pdir = getwd()
  if (verbose) {
    cat("Variables loading \n")
  }
  # setwd(path)
  Env <- stack()

  # Rasters loading
  files.null <- files
  if (is.null(files)) {
    files.null <- TRUE
  } else {
    files.null <- FALSE
  }
  resostack <- c(0, 0)
  extentstack <- extent(-Inf, Inf, -Inf, Inf)
  if (files.null) {
    n <- length(format)
  } else {
    n <- 1
  }
  for (j in 1:n) {
    if (files.null) {
      files <- list.files(path = path, pattern = paste0(".", format[j],
                                                        "$"))
    }
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        if (!is.null(path)) {
          file <- paste0(path, "/", files[[i]])
        } else {
          file <- files[i]
        }
        Raster <- raster(file)
        # Extent and resolution check
        reso <- res(Raster)
        extent <- extent(Raster)
        resostack[1] <- max(reso[1], resostack[1])
        resostack[2] <- max(reso[2], resostack[2])
        extentstack@xmin <- max(extentstack@xmin, extent@xmin)
        extentstack@xmax <- min(extentstack@xmax, extent@xmax)
        extentstack@ymin <- max(extentstack@ymin, extent@ymin)
        extentstack@ymax <- min(extentstack@ymax, extent@ymax)
        if (GUI) {
          incProgress(1/(length(files) * 3), detail = paste(i,
                                                            "loaded"))
        }
      }
    }
  }

  if (verbose) {
    cat("Variables treatment \n")
  }

  if (verbose) {
    cat("   resolution and extent adaptation...")
  }
  for (j in 1:n) {
    if (files.null) {
      files <- list.files(path = path, pattern = paste0(".", format[j],
                                                        "$"))
    }
    if (length(files) > 0) {
      for (i in seq_len(length(files))) {
        if (!is.null(path)) {
          file <- paste0(path, "/", files[[i]])
        } else {
          file <- files[i]
        }
        Raster <- raster(file)
        Raster <- readAll(Raster)
        Raster[Raster[] <= -900] <- NA
        names(Raster) <- as.character(strsplit(files[i], format[j]))
        if (names(Raster) %in% categorical) {
          Raster <- raster::as.factor(Raster)
          row.names(Raster@data@attributes[[1]]) <- Raster@data@attributes[[1]]$ID
          fun <- max
        } else {
          fun <- mean
        }
        if (round(res(Raster)[1], digits = 6) != round(resostack[1],
                                                       digits = 6) || round(res(Raster)[2], digits = 6) != round(resostack[2],
                                                                                                                 digits = 6)) {
          cat("\n", names(Raster), "is aggregated to the coarsest resolution \n")
          Raster <- aggregate(Raster, fact = c((resostack[1]/res(Raster)[1]),
                                               (resostack[2]/res(Raster)[2])), fun = fun)
        }
        Raster <- crop(Raster, extentstack)
        if (i > 1 && (any(res(Raster) != res(Env)) || extent(Raster) !=
                      extentstack)) {
          Raster <- resample(Raster, Env[[1]])
        }
        Env <- stack(Env, Raster, quick = TRUE)
        if (GUI) {
          incProgress((1/(length(files) * 3)), detail = paste(i,
                                                              "treated"))
        }
      }
    }
  }
  if (verbose) {
    cat("   done... \n")
  }


  # Normalizing variable
  if (Norm) {
    if (verbose) {
      cat("   normalizing continuous variable \n\n")
    }
    for (i in seq_len(length(Env@layers))) {
      # For not categorical variable
      if (!Env[[i]]@data@isfactor) {
        Env[[i]] <- Env[[i]]/Env[[i]]@data@max
      }
      if (GUI) {
        incProgress((1/length(Env@layers)/3), detail = paste(i, "normalized"))
      }
    }
  }

  # setwd(pdir)

  # Temporary files
  if (tmp) {
    path <- get("tmpdir", envir = .PkgEnv)

    if (!(file.path(path, ".rasters") %in% list.dirs(path)))
      (dir.create(paste0(path, "/.rasters")))

    if(!is.null(folder_tmp)){
      if (!(file.path(path, ".rasters", folder_tmp ) %in% list.dirs(paste0(path, "/.rasters"))))
        (dir.create(paste0(path, "/.rasters/",folder_tmp)))
      for (i in seq_len(length(Env@layers))) {
        Env[[i]] <- writeRaster(Env[[i]], paste0(path, "/.rasters/", folder_tmp, "/",
                                                 names(Env[[i]])), overwrite = TRUE)
      }
    }
    else{
      for (i in seq_len(length(Env@layers))) {
        Env[[i]] <- writeRaster(Env[[i]], paste0(path, "/.rasters/",
                                                 names(Env[[i]])), overwrite = TRUE)
      }
    }
  }

  return(Env)
}
