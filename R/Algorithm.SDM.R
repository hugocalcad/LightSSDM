#' @include SDM.R
#' @import methods
NULL

#' An S4 class to represent an SDM based on a single algorithm
#' \strong{(Esp)} Una clase S4 para representas un SDM basado para un algoritmo
#'
#' This is an S4 class to represent an SDM based on a single algorithm (including
#' generalized linear model, general additive model, multivariate adpative
#' splines, generalized boosted regression model, classification tree analysis,
#' random forest, maximum entropy, artificial neural network, and support vector
#' machines). This S4 class is obtained with \code{\link{modelling}}.
#' \strong{(Esp)} Es una clase S4 que representa un SDM basado sobre un algoritmo(Incluye
#' modelos lineales generalizados, modelos aditivos generalizados, reglas multivariadas
#' de regresión adaptativa, modelo aumentado de regresión generalizada, análisis de
#' clasificación de árboles,bosques aleatorios, entropia máxima, redes neuronales
#' artificiales , y máquinas vectoriales de apoyo. Esta clase S4 es obtenida
#' con \code{\link{modelling}}.
#'
#' @slot name character. Name of the SDM (by default Species.SDM).
#' \strong{(Esp)} Nombre de el SDM (por defecto Species.SDM)
#' @slot projection raster. Habitat suitability map produced by the SDM.
#' \strong{(Esp)} Mapa de idoneidad del hábitat producido por el SDM.
#' @slot binary raster. Presence/Absence binary map produced by the SDM.
#' \strong{(Esp)} Mapa bianrio de Presencia/Ausencia producido por el SDM
#' @slot evaluation data frame. Evaluation of the SDM (available metrics include
#'  AUC, Kappa, sensitivity, specificity and proportion of correctly predicted
#'  occurrences) and identification of the optimal threshold to convert the
#'  habitat suitability map into a binary presence/absence map.
#'  \strong{(Esp)} Evaluación de el SDM (las métricas disponibles incluyen
#'  AUC, Kappa, senibilidad, especificidad y proporción de ocurrencias
#'  correctamente predichos) e identificación del umbral óptimo para
#'  convertir el mapa de aptitud del hábitat en un mapa binario
#'  de presencia/ausencia.
#' @slot variable.importance data frame. Relative importance of
#'  each variable in the SDM.
#'  \strong{(Esp)} Importancia relativa de cada variabla en el SDM
#' @slot data data frame. Data used to build the SDM.
#' \strong{(Esp)} Datos usados para construir el SDM.
#' @slot parameters data frame. Parameters used to build the SDM.
#' \strong{(Esp)} Parámetros usados para construir el SDM.
#'
#' @seealso \linkS4class{Ensemble.SDM} an S4 class for ensemble SDMs,
#'  and \linkS4class{Stacked.SDM} an S4 class for SSDMs.
#'  \strong{(Esp)} \linkS4class{Ensemble.SDM} una clase S4 para ensamblar SDMs,
#'  y \linkS4class{Stacked.SDM} una clase S4 para SSDMs.
#'
#' @export
setClass('Algorithm.SDM',
         contains = 'SDM')

# Class generator
Algorithm.SDM <- function(algorithm = 'Algorithm',
                          name = character(),
                          projection = raster(),
                          binary = raster(),
                          evaluation = data.frame(),
                          variable.importance = data.frame(),
                          data = data.frame(),
                          parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  object.class <- paste0(algorithm,'.SDM')
  return(new(object.class,
             name = name,
             binary = binary,
             projection = projection,
             evaluation = evaluation,
             variable.importance = variable.importance,
             data = data,
             parameters = parameters))
}

setClass('GLM.SDM',
         contains = 'Algorithm.SDM')

setClass('BAM.SDM',
         contains = 'Algorithm.SDM')

setClass('MARS.SDM',
         contains = 'Algorithm.SDM')

setClass('CTA.SDM',
         contains = 'Algorithm.SDM')

setClass('GBM.SDM',
         contains = 'Algorithm.SDM')

setClass('RF.SDM',
         contains = 'Algorithm.SDM')

setClass('MAXNET.SDM',
         contains = 'Algorithm.SDM')

setClass('ANN.SDM',
         contains = 'Algorithm.SDM')

setClass('KSVM.SDM',
         contains = 'Algorithm.SDM')
