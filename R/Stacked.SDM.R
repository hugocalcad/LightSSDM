#' @include Ensemble.SDM.R
#' @importFrom raster raster stack
NULL

#'An S4 class to represent SSDMs
#'\strong{(Esp)} Una clase S4 que represta SSDMs
#'
#'This is an S4 class to represent SSDMs that assembles multiple algorithms
#'(including generalized linear model, general additive model, multivariate
#'adaptive splines, generalized boosted regression model, classification tree
#'analysis, random forest, maximum entropy, artificial neural network, and
#'support vector machines) built for multiple species. It is obtained with
#'\code{\link{stack_modelling}} or \code{\link{stacking}}.
#'\strong{(Esp)} Es una clase S4 para representar SSDM que ensambla múltiples
#'algoritmos (incluido el modelo lineal generalizado, el modelo aditivo general,
#'las splines adaptativas multivariadas, el modelo de regresión reforzada generalizada,
#'el análisis del árbol de clasificación, el bosque aleatorio, la entropía máxima,
#'la red neuronal artificial y las máquinas de vectores de soporte) Construido para
#'múltiples especies. Se obtiene con \code{\link{stack_modelling}} o \code{\link{stacking}}.
#'
#'@slot name character. Name of the SSDM (by default 'Species.SSDM').
#'\strong{(Esp)} Nombre del SSDM(por defecto 'Species.SDM')
#'@slot diversity.map raster. Local species richness map produced by the SSDM.
#'\strong{(Esp)} Mapa local de riqueza de especies producido por el SSDM.
#'@slot endemism.map raster. Endemism map produced by the SSDM (see Crisp et al
#'  (2011) in references).
#'\strong{(Esp)} Mapa de endemismo producido por el SSDM (ver Crisp et al (2011) en referencias)
#'@slot uncertainty raster. Between-algorithm variance map.
#'\strong{(Esp)} Mapa de varianza entre algoritmos
#'@slot evaluation data frame. Evaluation of the SSDM (AUC, Kappa, omission
#'  rate, sensitivity, specificity, proportion of correctly predicted occurrences).
#'\strong{(Esp)} Evaluacion del SSDM (AUC, Kappa, tasa de omisión, sesitividad, especificidad,
#'  proporción de ocurrencias correctamente predichas)
#'@slot variable.importance data frame. Relative importance of each variable in the SSDM.
#'\strong{(Esp)} Importancia relativa de cada variable en el SSDM.
#'@slot algorithm.correlation data frame. Between-algorithm correlation matrix.
#'\strong{(Esp)} Matriz de correlación entre algoritmos.
#'@slot enms list. List of ensemble SDMs used in the SSDM.
#'\strong{(Esp)} Lista de SDMs ensamblados en el SSDM.
#'@slot parameters data frame. Parameters used to build the SSDM.
#'\strong{(Esp)} Parametros usados para construir el SSDM.
#'@slot algorithm.evaluation data frame. Evaluation of the algorithm averaging
#'  the metrics of all SDMs (AUC, Kappa, omission rate, sensitivity,
#'  specificity, proportion of correctly predicted occurrences).
#'\strong{(Esp)} Evaluación del algoritmo que promedia las métricas de todos los SDM
#'  (AUC, Kappa, tasa de omisión, sensibilidad, especificidad, proporción de ocurrencias
#'  predichas correctamente).
#'
#'@seealso \linkS4class{Ensemble.SDM} an S4 class to represent ensemble SDMs,
#'  and \linkS4class{Algorithm.SDM} an S4 class to represent SDMs.
#'\strong{(Esp)} Un objeto de clase S4 que representa los SDMs ensamblados, y una clase
#'  S4 \linkS4class{Algorithm.SDM} que representa los SDMs.
#'
#'@references M. D. Crisp, S. Laffan, H. P. Linder & A. Monro (2001)
#'  "Endemism in the Australian flora"  \emph{Journal of Biogeography}
#'  28:183-198
#'  \url{http://biology-assets.anu.edu.au/hosted_sites/Crisp/pdfs/Crisp2001_endemism.pdf}
#'
#'
#'
#'@export
setClass('Stacked.SDM',
         representation(name = 'character',
                        diversity.map = 'Raster',
                        endemism.map = 'Raster',
                        uncertainty = 'Raster',
                        evaluation = 'data.frame',
                        variable.importance = 'data.frame',
                        algorithm.correlation = 'data.frame',
                        algorithm.evaluation = 'data.frame',
                        enms = 'list',
                        parameters = 'data.frame'),
         prototype(name = character(),
                   diversity.map = raster(),
                   endemism.map = raster(),
                   uncertainty = raster(),
                   evaluation = data.frame(),
                   variable.importance = data.frame(),
                   algorithm.correlation = data.frame(),
                   algorithm.evaluation = data.frame(),
                   enms = list(),
                   parameters = data.frame()))

# Class Generator
Stacked.SDM <- function(name = character(),
                        diversity.map = raster(),
                        endemism.map = raster(),
                        uncertainty = raster(),
                        evaluation = data.frame(),
                        variable.importance = data.frame(),
                        algorithm.correlation = data.frame(),
                        algorithm.evaluation = data.frame(),
                        enms = list(),
                        parameters = data.frame(matrix(nrow = 1, ncol = 0))) {
  return(new('Stacked.SDM',
             name = name,
             diversity.map = diversity.map,
             endemism.map = endemism.map,
             evaluation = evaluation,
             variable.importance = variable.importance,
             uncertainty = uncertainty,
             algorithm.correlation = algorithm.correlation,
             algorithm.evaluation = algorithm.evaluation,
             enms = enms,
             parameters = parameters))}

