#'@include ensemble_modelling.R stacking.R checkargs.R
NULL

#'Update a previous SSDM
#'\strong{(Esp)} Actualiza un SSDM previo
#'
#'Update a previous SSDM with new occurrence data. The function takes as inputs
#'updated or new occurrence data from one species, previous environmental
#'variables, and an S4 \linkS4class{Stacked.SDM} class object containing a
#'previously built SSDM.
#'\strong{(Esp)} Actualizar un SSDM anterior con nuevos datos de ocurrencia.
#'La función toma como entradas actualizadas o nuevos datos de ocurrencia de
#'una especie, variables ambientales anteriores y un objeto de clase S4
#'\linkS4class{Stacked.SDM} que contiene un SSDM creado anteriormente.
#'
#'@param object Stacked.SDM. The previously built SSDM.
#'\strong{(Esp)} El SSDM construido previamente.
#'@param Occurrences data frame. New or updated occurrence table (can be
#'  processed first by \code{\link{load_occ}}).
#'\strong{(Esp)} Tabla de ocurrencias nnueva o actualizada (puede ser procesada
#'primero por \code{\link{load_occ}}).
#'@param Env raster object. Environment raster object (can be processed first by
#'  \code{\link{load_var}}).
#'\strong{(Esp)} Objeto raster ambiental (puede sr procesado primero por \code{\link{load_var}})
#'@param Xcol character. Name of the column in the occurrence table containing
#'  Latitude or X coordinates.
#'\strong{(Esp)} Nombre de la columna en la tabla de ocurrencias que contiene
#' la latitud o la coordenada X.
#'@param Ycol character. Name of the column in the occurrence table containing
#'  Longitude or Y coordinates.
#'\strong{(Esp)} Nombre de la columna en la tabla de ocurrencias que contiene
#' la longitud o la coordenada Y.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  whether a line is a presence or an absence. A value of 1 is presence and
#'  value of 0 is absence. If NULL presence-only dataset is assumed.
#'\strong{(Esp)} Nombre de la columna en la tabla de ocurrencias especificando
#'  si la linea es una presencia o ausencia. Un valor de 1 es presencia y un valor
#'  de 0 es ausencia. Si es NULL se asume solo datos de presencia.
#'@param Spname character. Name of the new or updated species.
#'\strong{(Esp)} Nombre de la especie nuevaa o actualizada.
#'@param name character. Optional name given to the final SSDM produced, by
#'  default it's the name of the previous SSDM.
#'\strong{(Esp)} Nombre opcional dado al SDM final producido (por defecto es el
#'  nombre del SSDM previo)
#'@param save logical. If set to true, the model is automatically saved.
#'\strong{(Esp)} Si es verdadero (\code{TrUE}) el modelo es automaticamente guardado.
#'@param path character. Name of the path to the directory to contain the saved SSDM.
#'\strong{(Esp)} Nombre de la ruta del directorio que contiene el SSDM guardado.
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see
#'  \code{\link[SDMTools]{optim.thresh}}).
#'\strong{(Esp)} Un valor numérico simple  representado el numero de igual intervalo
#'  del umbral entre 0 y 1 (ver \code{\link[SDMTools]{optim.thresh}}).
#'@param tmp logical. If set to true, the habitat suitability map of each
#'  algorithm is saved in a temporary file to release memory. But beware: if you
#'  close R, temporary files will be deleted To avoid any loss you can save
#'  your model with \code{\link{save.model}}.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}) el mapa de idoneidad de hábitat de cada
#'  algoritmo es guardado en un archivo temporal para liberar memoria. Pero hay que
#'  tener cuidado: Si se cierra R, los archivos temporales seran borrados para evitar cualquier
#'  perdida puedes guardar tu modelo con \code{\link{save.model}}.
#'@param verbose logical. If set to true, allows the function to print text in the console.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}), permite la función de imprimir texto en la consola.
#'@param GUI logical. Don't take that argument into account (parameter for the user interface).
#'\strong{(Esp)} No tomar este arguento en cuenta(Parametro para la interfaz de usuario)
#'@param ... additional parameters for the algorithm modelling function (see details below).
#'\strong{(Esp)} parametros adicionales la funcion modelling de algoritmos (ver detalles más abajo)
#'
#'@return an S4 \linkS4class{Stacked.SDM} class object viewable with the
#'  \code{\link{plot.model}} function.
#'\strong{(Esp)} Un objeto de clase S4 \linkS4class{Stacked.SDM} visible con la función
#'  \code{\link{plot.model}}.
#'
#'@seealso \code{\link{stack_modelling}} to build SSDMs.
#'\strong{(Esp)} Ver también \code{\link{stack_modelling}} para construir SSDMs.
#'
#' @examples
#' \dontrun{
#' update(stack, Occurrences, Env, Spname = 'NewSpecie')
#' }
#'
#'@export
setMethod('update', 'Stacked.SDM',
          function(object,
                   # Modelling data input
                   Occurrences, Env,
                   # Occurrences reading
                   Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL, Spname = NULL,
                   # Model creation
                   name = stack@name, save = FALSE, path = getwd(), thresh = 1001, tmp = FALSE,
                   # Verbose
                   verbose = TRUE, GUI = FALSE,
                   # Modelling parameters
                   ...) {
            # Check arguments
            .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, Spname = Spname, name = name,
                       save = save, path = path, thresh = thresh, tmp = tmp, verbose = verbose,
                       GUI = GUI)

            stack <- object
            # New ENM creation
            if (verbose) {
              cat("New species ensemble distribution model creation...\n")
            }
            if (stack@parameters$PA == "default") {
              PA <- NULL
            } else {
              PA <- list(nb = strsplit(stack@parameters$PA, ".", fixed = TRUE)[[1]][1],
                         strat = strsplit(stack@parameters$PA, ".", fixed = TRUE)[[1]][2])
            }
            if (!is.null(Spname)) {
              enm.name <- Spname
            } else {
              enm.name <- "new_Specie"
            }
            ENM <- ensemble_modelling(strsplit(stack@parameters$algorithms, ".", fixed = TRUE)[[1]][-1],
                                      Occurrences, Env, Xcol, Ycol, Pcol, rep = as.numeric(stack@parameters$rep),
                                      enm.name, save = FALSE, path = getwd(), PA, cv = stack@parameters$cv,
                                      cv.param = as.numeric(strsplit(stack@parameters$cv.param, "|", fixed = TRUE)[[1]][-1]),
                                      thresh = thresh, metric = stack@parameters$metric, axes.metric = stack@parameters$axes.metric,
                                      uncertainity = stack@uncertainity@data@haveminmax, tmp = tmp,
                                      ensemble.metric = strsplit(stack@parameters$ensemble.metric,
                                                                 ".", fixed = TRUE)[[1]][-1],
                                      ensemble.thresh = as.numeric(strsplit(stack@parameters$ensemble.thresh,
                                                                            "|", fixed = TRUE)[[1]][-1]),
                                      weight = as.logical(stack@parameters$weight),
                                      ...)
            if (verbose) {
              cat("   done.\n")
            }

            # Test for new
            if (verbose) {
              cat("Check if the species already exist...\n")
            }
            if (!is.null(Spname)) {
              i <- which(names(stack@enms) == paste0(Spname, ".Ensemble.SDM"))
              if (verbose) {
                cat(Spname, "replacement\n")
              }
              if (length(i) > 0) {
                stack@enms[[i]] <- NULL
              } else {
                stack@parameters$sp.nb.origin <- stack@parameters$sp.nb.origin +
                  1
              }
            } else {
              stack@parameters$sp.nb.origin <- stack@parameters$sp.nb.origin + 1
            }
            if (verbose) {
              cat("   done.\n")
            }

            # New stacking
            if (verbose) {
              cat("New stacking...\n")
            }
            enms <- list()
            for (i in seq_len(stack@enms)) {
              enms[[i]] <- stack@enms[[i]]
            }
            enms["method"] <- stack@parameters$method
            enms["endemism"] <- strsplit(stack@parameters$endemism, "|", fixed = "T")[[1]]
            enms["rep.B"] <- stack@parameters$rep.B
            newstack <- do.call(stacking, enms)
            if (verbose) {
              cat("   done.\n")
            }

            if (!is.null(stack)) {
              # Paremeters
              newstack@parameters$sp.nb.origin <- stack@parameters$sp.nb.origin

              # Saving
              if (save) {
                if (verbose) {
                  cat("Saving...\n")
                }
                if (!is.null(name)) {
                  save.stack(newstack, name = name, path = path)
                } else {
                  save.stack(newstack, path = path)
                }
                if (verbose) {
                  cat("   done.\n")
                }
              }
            }

            return(newstack)
          })
