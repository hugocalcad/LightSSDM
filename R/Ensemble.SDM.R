#' @include Algorithm.SDM.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent an ensemble SDM
#'\strong{(Esp)} Una calse S4 que representa un SDM ensamblado
#'
#'This is an S4 class to represent an ensemble SDM from multiple algorithms
#'(including generalized linear model, general additive model, multivariate
#'adaptive splines, generalized boosted regression model, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines). This S4 class is returned by
#'\code{\link{ensemble_modelling}} or \code{\link{ensemble}}.
#'\strong{(Esp)} Esta es una clase S4 que representa un SDm ensamblado de multiples algoritmos
#'(incluye modelos lineales generalizados, modelos aditivos generalizados, reglas multivariadas
#' de regresión adaptativa, modelo aumentado de regresión generalizada, análisis de
#' clasificación de árboles,bosques aleatorios, entropia máxima, redes neuronales
#' artificiales , y máquinas vectoriales de apoyo). Esta clase es devuelta por
#' \code{\link{ensemble_modelling}} ó \code{\link{ensemble}}.
#'
#'@slot uncertainty raster. Between-algorithm variance map.
#'\strong{(Esp)} Mapa de varianza entre algoritmos.
#'@slot algorithm.correlation data frame. Between-algorithm correlation matrix.
#'\strong{(Esp)} Matriz de correlación entre algoritmos.
#'@slot algorithm.evaluation data frame. Evaluation of the ensemble SDM (available
#'  metrics include AUC, Kappa, sensitivity, specificity and proportion of
#'  correctly predicted occurrences) and identification of the optimal threshold
#'  to convert the habitat suitability map into a binary presence/absence map.
#'\strong{(Esp)} Evaluación de el SDM ensamblado (métricas disponibles incluido
#'   AUC, Kappa, sensitividad, especifidad y proporción de predicción de ocurrencias
#'  correcta) e identificación del umbral optimo para convertir el mapa de habitabilidad
#'  en un mapa binario de presencia/ausencia.
#'
#'@seealso \linkS4class{Algorithm.SDM} an S4 class to represent an SDM based on
#'  a single algorithm, and \linkS4class{Stacked.SDM} an S4 class for SSDMs.
#'\strong{(Esp)} \linkS4class{Algorithm.SDM} un clase S4 para representas un SDM
#'Basado en un solo algoritmo, y \linkS4class{Stacked.SDM} una clase S4 para SSDMs.
#'
#'@export
setClass('Ensemble.SDM',
         contains = 'SDM',
         representation(uncertainty = 'Raster',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame'),
         prototype(uncertainty = raster(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame()))

Ensemble.SDM <- function(name = character(),
                         projection = raster(),
                         binary = raster(),
                         evaluation = data.frame(),
                         variable.importance = data.frame(),
                         data = data.frame(),
                         uncertainty = raster(),
                         algorithm.correlation = data.frame(),
                         algorithm.evaluation = data.frame(),
                         parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Ensemble.SDM',
             name = name,
             projection = projection,
             binary = binary,
             evaluation = evaluation,
             variable.importance = variable.importance,
             data = data,
             uncertainty = uncertainty,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             parameters = parameters))}
