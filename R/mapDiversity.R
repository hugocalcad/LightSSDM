#' @include Ensemble.SDM.R checkargs.R
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPoints bbox
#' @importFrom raster raster stack reclassify mask calc overlay values rasterize rasterToPoints values<-
#' @importFrom stats lm
NULL

#' Map Diversity
#' \strong{(Esp)} Mapa de diversidad
#'
#' Methods for Stacked.SDM or SSDM to map diversity and communities composition.
#' \strong{(Esp)} Métodos para Stacked.SDM ó SSDM para mapas de diversidad y composición decomunidades.
#'
#'@param obj Stacked.SDM. SSDM to map diversity with.
#'  \strong{(Esp)} SSDM con el que se va a mapear la diversidad
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'  \strong{(Esp)} define el método usado para crear la distribución local de especies.
#'  (ver detalles abajo)
#'@param rep.B integer. If the method used to create the local species richness
#'  is the random Bernoulli (\strong{Bernoulli}), rep.B parameter defines the number of
#'  repetitions used to create binary maps for each species.
#'  \strong{(Esp)} Si el método usado para crear la distribución local de especies es el Bernoulli
#'  aleatorio (\strong{Bernoulli}), se define como el númeor de repeticiones para crear el mapa
#'  binario por cada especie.
#'@param Env raster object. Stacked raster object of environmental variables
#'  (can be processed first by \code{\link{load_var}}). Needed only for stacking
#'  method using probability ranking from richness (\strong{PRR}).
#'  \strong{(Esp)} Pila de rasters de las variables ambientales (que puede ser procesado primero
#'  por \code{\link{load_var}})). Solo se necesita para el método de apilamiento mediante la
#'  clasificación de probabilidad de la riqueza (PRR)
#'@param verbose logical. If set to true, allows the function to print text in
#'  the console.
#'  \strong{(Esp)} Si es verdadero (\code{TRUE}) permite imprimir texto en la consola.
#'@param ... other arguments pass to the method.
#'  \strong{(Esp)} otros argumentos
#'
#'@return a list with a diversity map and eventually ESDMs for stacking method
#'  using probability ranking from richness (\strong{PPR}).
#'  \strong{(Esp)} Una lista con un mpara de diversidad y eventuales ESDMs para el método de apilamiento mediante la
#'  clasificación de probabilidad de la riqueza (PRR)
#'
#'@details \strong{Methods:} Choice of the method used to compute the local
#'  species richness map (see Calabrese et al. (2014) and D'Amen et al (2015) for
#'  more informations, see reference below): \strong{(Esp)} Elige el método usado para procesar
#'  el mapa de riqueza de especies (ver Calabrese et al. (2014) y D'Amen el al (2015) para
#'  más información ver referencias abajo): \describe{\item{pSSDM}{sum
#'  probabilities of habitat suitability maps. \strong{(Esp)} suma de probabilidades de los
#'  mapas de idoneidad del hábitat}\item{Bernoulli}{draw repeatedly
#'  from a Bernoulli distribution. \strong{(Esp)} dibujar repetidamente una distribución de Bernoulli}
#'  \item{bSSDM}{sum the binary map obtained with the thresholding (depending on the metric of the
#'  ESDM). \strong{(Esp)} suma el mapa binario obtenido con el umbral (según la métrica del ESDM).}
#'  \item{MaximumLikelihood}{adjust species richness of the model by
#'  linear regression. \strong{(Esp)} Ajustar la riqueza de especies del modelo mediante regresión lineal.}
#'  \item{PRR.MEM}{model richness with a macroecological model (MEM) and adjust each ESDM binary map by
#'  ranking habitat suitability and keeping as much as predicted richness of the MEM. \strong{(Esp)}
#'  Modela la riqueza con un modelo macroecológico (MEM) y ajuste cada mapa binario de ESDM clasificando
#'  la aptitud del hábitat y manteniendo la riqueza predicha del MEM.}\item{PRR.pSSDM}{model
#'  richness with a pSSDM and adjust each ESDM binary map by ranking habitat
#'  suitability and keeping as much as predicted richness of the pSSDM. \strong{(Esp)} Modela la riqueza con un
#'  pSSDM y ajuste cada mapa binario de ESDM clasificando la idoneidad del hábitat y manteniendo la riqueza predicha del pSSDM}}
#'
#'@examples
#'
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' # SSDM building
#' SSDM <- stack_modelling(c('CTA', 'KSVM'), Occurrences, Env, rep = 1,
#'                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                        Spcol = 'SPECIES')
#'
#' # Diversity mapping
#' mapDiversity(SSDM, mathod = 'pSSDM')
#'
#' }
#'
#'@seealso \code{\link{stacking}} to build SSDMs.
#' \strong{(Esp)} ver \code{\link{stacking}} para construir SSDMs.
#'
#'@references M. D'Amen, A. Dubuis, R. F. Fernandes, J. Pottier, L. Pelissier, &
#'  A Guisan (2015) "Using species richness and functional traits prediction to
#'  constrain assemblage predicitions from stacked species distribution models"
#'  \emph{Journal of Biogeography} 42(7):1255-1266
#'  \url{http://doc.rero.ch/record/235561/files/pel_usr.pdf}
#'
#'  J.M. Calabrese, G. Certain, C.  Kraan, & C.F. Dormann (2014) "Stacking
#'  species distribution  models  and  adjusting  bias  by linking them to
#'  macroecological models." \emph{Global Ecology and Biogeography} 23:99-112
#'  \url{http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/calabrese2013globalecolbiogeogr.pdf}
#'
#' @name mapDiversity
NULL

#' @rdname mapDiversity
#' @export
setGeneric('mapDiversity', function(obj, ...) {return(standardGeneric('mapDiversity'))})

#' @rdname mapDiversity
#' @export
setMethod("mapDiversity", "Stacked.SDM", function(obj, method, rep.B = 1000,
                                                  verbose = TRUE, Env = NULL,
                                                  ...){
  # Check arguments
  .checkargs(stack = obj, method = method, rep.B = rep.B, verbose = verbose)

  enms <- NULL # Preparing enms slot for PPR methods
  diversity.map <- reclassify(obj@enms[[1]]@projection[[1]], c(-Inf,Inf, 0))

  # Useless datacheck to prevent bugs to remove after debugging
  for (i in seq_len(length(obj@enms))) {
    if (!inherits(obj@enms[[i]]@projection, "RasterLayer")) {
      if (verbose) {
        cat("Error", obj@enms[[i]]@name, "is not a raster but a",
            class(obj@enms[[i]]@projection)[1],
            ".\nIt will be removed for the stacking")
      }
      obj@enms[[i]] <- NULL
    }
  }

  if (method == "bSSDM") {
    # Threshold and sum (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by thresholding and then summing. \n")
    diversity.map <- sum(stack(lapply(obj@enms, function(x)
      reclassify(x@projection,
                 c(-Inf,x@evaluation$threshold,0,x@evaluation$threshold,Inf,1))
      )))
  }

  if (method == "pSSDM") {
    # Individual probabilities sum (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by summing individual probabilities. \n")
    diversity.map <- sum(stack(lapply(obj@enms, function(x) x@projection)))
  }

  if (method == "Bernoulli") {
    # Random Bernoulli distribution (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by drawing repeatedly from a Bernoulli distribution. \n")
    proba <- stack(lapply(obj@enms, function(x) x@projection))
    diversity.map <- calc(proba, fun = function(...) {
      x <- c(...)
      x[is.na(x)] <- 0
      return(rbinom(lengths(x), rep.B, x))
    }, forcefun = TRUE)
    diversity.map <- sum(diversity.map)/length(enms)/rep.B
  }

  if (method == "MaximumLikelihood") {
    # Maximum likelihood (Calabrese et al, 2014)
    if (verbose)
      cat("\n Local species richness computed by maximum likelihood adjustment. \n")
    diversity.map <- mapDiversity(obj, method = 'bSSDM',
                                  verbose = FALSE)$diversity.map
    Richness <- .richness(obj)
    SSDM_Richness <- values(mask(diversity.map, Richness))
    SSDM_Richness <- SSDM_Richness[-which(is.na(SSDM_Richness))]
    Richness <- values(Richness)
    Richness <- Richness[-which(is.na(Richness))]
    fit <- lm(Richness ~ SSDM_Richness)
    a <- fit$coefficients[1]
    b <- fit$coefficients[2]
    diversity.map <- a + b * diversity.map
  }

  if (method == "PRR.MEM") {
    # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
    if (verbose)
      cat("\n Local species richness computed by probability ranking from MEM. \n")
    diversity.map <- .MEM(obj,Env)@projection
    enms <- .PRR(obj, diversity.map)
  }

  if (method == "PRR.pSSDM") {
    # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
    if (verbose)
      cat("\n Local species richness computed by probability ranking from pSSDM. \n")
    diversity.map <- mapDiversity(obj, method = 'pSSDM',
                                  verbose = FALSE)$diversity.map
    enms <- .PRR(obj, diversity.map)
  }

  return(list(
    diversity.map = diversity.map,
    enms = enms
  ))
})

##### Internals ####

.richness <- function(obj){
  Richness <- reclassify(obj@enms[[1]]@projection, c(-Inf, Inf, 0))
  for (i in seq_len(length(obj@enms)))
    Richness <- Richness + rasterize(
      SpatialPoints(obj@enms[[i]]@data[1:2]),
      Richness, field = obj@enms[[i]]@data$Presence,
      background = 0)
  if (all(values(Richness) %in% c(0, 1, NA)))
    stop("Observed Richness is always equal to 1, modelled richness can't be adjusted !")
  return(Richness)
}

.MEM <- function(obj, Env){
  occ <- data.frame(rasterToPoints(.richness(obj), function(x) x > 0))
  maxOcc <- max(occ$layer) # Reucing occ for algorithms
  occ$layer <- occ$layer/max(maxOcc)
  algo <- unlist(
    strsplit(obj@enms[[1]]@parameters$algorithms,
             ".", fixed = TRUE))[-1]
  MEM <- ensemble_modelling(algorithms = algo,
    Occurrences = occ, Env = Env, Xcol = "x",
    Ycol = "y", Pcol = "layer", rep = obj@enms[[1]]@parameters$rep,
    name = "MEM", cv = obj@enms[[1]]@parameters$cv,
    cv.param = as.numeric(unlist(
      strsplit(obj@enms[[1]]@parameters$cv.param,
               "|", fixed = TRUE))[-1]),
    metric = obj@enms[[1]]@parameters$metric,
    axes.metric = obj@enms[[1]]@parameters$axes.metric,
    ensemble.metric = unlist(
      strsplit(obj@enms[[1]]@parameters$ensemble.metric,
               ".", fixed = TRUE))[-1],
    ensemble.thresh = as.numeric(unlist(
      strsplit(obj@enms[[1]]@parameters$ensemble.thresh,
               "|", fixed = TRUE))[-1]),
    uncertainty = FALSE,
    weight = as.logical(obj@enms[[1]]@parameters$weight),
    verbose = FALSE)
  MEM@projection <- MEM@projection*maxOcc
  return(MEM)
}

.PRR <- function(obj, Richness){
  # Readjust each enm binary map
  richnesses <- values(Richness)
  names(richnesses) <- seq_len(length(richnesses))
  richnesses <- as.list(richnesses)
  probabilities <- lapply(lapply(obj@enms, FUN = slot, name = "projection"),
                          values)
  probabilities <- lapply(probabilities, function(x) {
    names(x) <- rep(seq_len(length(probabilities[[1]])))
    return(x)
  })
  probabilities <- lapply(probabilities, `[`, names(probabilities[[1]]))
  probabilities <- apply(do.call(rbind, probabilities), 2, as.list)
  binaries <- lapply(lapply(obj@enms, FUN = slot, name = "binary"),
                     values)
  binaries <- lapply(binaries, function(x) {
    names(x) <- rep(seq_len(length(binaries[[1]])))
    return(x)
  })
  binaries <- lapply(binaries, `[`, names(binaries[[1]]))
  binaries <- apply(do.call(rbind, binaries), 2, as.list)
  binaries <- mapply(function(rich, probability, binary) {
    if (!is.na(rich)) {
      ord <- order(unlist(probability), decreasing = TRUE)
      binary <- unlist(binary)
      if (length(ord) <= rich) {
        binary[ord] <- 1
      } else {
        binary[ord[1:rich]] <- 1
        binary[ord[rich + seq_len(length(ord))]] <- 0
      }
      binary <- as.list(binary)
    }
    return(binary)
  }, rich = richnesses, probability = probabilities, binary = binaries,
  SIMPLIFY = FALSE)
  binaries <- lapply(binaries, `[`, names(binaries[[1]]))
  binaries <- apply(do.call(rbind, binaries), 2, as.list)
  binaries <- lapply(binaries, unlist)
  binaries <- lapply(binaries, unname)
  mapply(function(enm, binary) {
    values(enm@binary) <- binary
    return(enm)
  }, enm = obj@enms, binary = binaries, SIMPLIFY = FALSE)
}
