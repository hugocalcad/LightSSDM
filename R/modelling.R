#' @include Algorithm.SDM.R checkargs.R
#' @importFrom shiny incProgress
#' @importFrom raster stack writeRaster
NULL

#'Build an SDM using a single algorithm
#'\strong{(Esp)} Contruye un SDM usando un algoritmo simple
#'
#'This is a function to build an SDM with one algorithm for a single species.
#'The function takes as inputs an occurrence data frame made of presence/absence
#'or presence-only records and a raster object for data extraction and
#'projection. The function returns an S4 \linkS4class{Algorithm.SDM} class
#'object containing the habitat suitability map, the binary map and the
#'evaluation table.
#'\strong{(Esp)} Esta es una función para construir un SDM con un solo algoritmo para una sola especie.
#'La funcion como entrada una tabla de ocurrencias hecha de presencias/ausencias o
#'precias solamente y un objeto ráster para extracción de datos y projection. La función
#'regresa  un objeto de clase S4 \linkS4class{Algorithm.SDM} que contiene el mapa de
#'habitabilidad, el mapa binario y una tabla de evaluación.
#'
#'@param algorithm character. Choice of the algorithm to be run (see details below).
#'  \strong{(Esp)} Escoge el algoritmo a ser ejecutado (ver detalles más abajo)
#'@param Occurrences data frame. Occurrence table (can be processed first by \code{\link{load_occ}}).
#' \strong{(Esp)} Tabla de ocurrencias (pueden ser procesadas primero por \code{\link{load_occ}}).
#'@param Env raster object. Raster object of environmental variable (can be
#'  processed first by \code{\link{load_var}}).
#'  #'\strong{(Esp)} Objeto RasterStack de las variables ambientales(pueden ser
#'  procesados primero por \code{\link{load_var}})
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
#'@param name character. Optional name given to the final SDM produced (by default 'Algorithm.SDM').
#'  \strong{(Esp)} Nombre opcional dado al SDM final producido (por defecto 'Algorithm.SDM')
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only dataset. If PA is NULL, recommended PA selection
#'  strategy is used depending on the algorithms (see details below).
#'\strong{(Esp)} definieniendo la estrategia de selección de las pseudo-ausencias usadas
#'  en caso de datos con solo presencias. Si PA es nulo (NULL), estrategias recomendadas
#'  son usadas dependiendo del algoritmo (ver detalles más abajo)
#'@param cv character. Method of cross-validation used to evaluate the SDM (see details below).
#'\strong{(Esp)} Metodo de validación-cruzada usada para evaluar el SDM (ver detalles más abajo)
#'@param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the SDM (see details below).
#'\strong{(Esp)} Parametros asociados al método de validación- cruzada usada
#'  para evaluar el SDM (ver detalles más abajo).
#'@param select logical. If set to true, models are evaluated before being
#'  projected, and not kept if they don't meet selection criteria (see details below).
#'\strong{(Esp)}. Si es verdadero (\code{TRUE}), el modelo es evaluado antes de ser proyectado,
#'  y no se conservan si no cumplen los criterios de selección (ver detalles más abajo).
#'@param select.metric character. Metric(s) used to pre-select SDMs that reach a
#'  sufficient quality (see details below).
#'\strong{(Esp)} Métrica(s) usada(s) para preselecciona SDMs que lleguen a tener suficiente
#'  calidad (Ver detalles más abajo)
#'@param select.thresh numeric. Threshold(s) associated with the metric(s) used
#'  to compute the selection.
#'\strong{(Esp)} Umbral(es) asociado a la(s) métrica(s) utilizada(s) para calcular la selección.
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see \code{\link[SDMTools]{optim.thresh}}).
#'\strong{(Esp)} Un valor numérico simple  representado el numero de igual intervalo
#'  del umbral entre 0 y 1 (ver \code{\link[SDMTools]{optim.thresh}}).
#'@param metric character. Metric used to compute the binary map threshold (see details below).
#'\strong{(Esp)} Métrica usada para procesar el umbral del mapa binario(Ver más detalles abajo)
#'@param axes.metric Metric used to evaluate variable relative importance (see details below).
#'\strong{(Esp)} Métrica usada para evaluar la variable relativa de importancia (ver detalles más abajo)
#'@param verbose logical. If set to true, allows the function to print text in the console.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}), permite la función de imprimir texto en la consola.
#'@param GUI logical. Don't take that argument into account (parameter for the user interface).
#'\strong{(Esp)} No tomar este arguento en cuenta(Parametro para la interfaz de usuario)
#'@param ... additional parameters for the algorithm modelling function (see details below).
#'\strong{(Esp)} parametros adicionales la funcion modelling de algoritmos (ver detalles más abajo)
#'
#'@return an S4 \linkS4class{Algorithm.SDM} Class object viewable with the \code{\link{plot.model}} method.
#'\strong{(Esp)} un objeto de clase S4 \linkS4class{Algorithm.SDM} visble con el método \code{\link{plot.model}}
#'
#'@details \describe{ \item{algorithm}{'all' calls all the following
#'  algorithms. Algorithms include Generalized linear model (\strong{GLM}),
#'  Generalized additive model for large dataSet splines (\strong{BAM}),
#'  Multivariate adaptive regression  (\strong{MARS}), Boosted regressions model
#'  (\strong{GBM}), Classification tree analysis (\strong{CTA}), Random forest
#'  (\strong{RF}), Maximum entropy (\strong{MAXNET}), Artificial neural network
#'  (\strong{ANN}), and Support vector machines (\strong{KSVM}). Each algorithm
#'  has its own parameters settable with the \strong{...} (see each algorithm
#'  section below to set their parameters).
#'  \strong{(Esp)} 'all' incluye todos los algoritmos. Algoritmos incluidos modelos lineales
#'  generalizados (\strong{GLM}), modelos aditivos generalizados para datos grandes(\strong{BAM}),
#'  reglas multivariadas de regresión adaptativa (\strong{MARS}), modelo aumentado de regresión
#'  generalizada (\strong{GBM}), análisis de clasificación de árboles(\strong{CTA}),bosques aleatorios
#'  (\strong{RF}), entropia máxima (\strong{MAXNET}), redes neuronales artificiales (\strong{ANN}),
#'  y máquinas vectoriales de apoyo (\strong{KSVM})} \item{"PA"}{list with two values:
#'  \strong{nb} number of pseudo-absences selected, and \strong{strat} strategy
#'  used to select pseudo-absences: either random selection or disk selection.
#'  We set default recommendation from Barbet-Massin et al. (2012) (see
#'  reference).
#'  Lista con dos valores: \strong{nb} número de pseudo-ausencias seleccionados, y \strong{strat}
#'  estrategia usada para seeccionar pseudo- ausencias: selección aleatoria o selección de disco.
#'  Se colocó parametros recomendadospor Barbet-Massin et al. (2012) (ver referencias)}
#'  \item{cv}{\strong{Cross-validation} method used to split the
#'  occurrence dataset used for evaluation: \strong{holdout} data are
#'  partitioned into a training set and an evaluation set using a fraction
#'  (\emph{cv.param[1]}) and the operation can be repeated (\emph{cv.param[2]})
#'  times, \strong{k-fold} data are partitioned into k (\emph{cv.param[1]})
#'  folds being k-1 times in the training set and once the evaluation set and
#'  the operation can be repeated (\emph{cv.param[2]}) times, \strong{LOO}
#'  (Leave One Out) each point is successively taken as evaluation data.
#'  \strong{(Esp)} \strong{Validación-cruzada} método utilizado para dividir el conjunto de datos de las
#'  ocurrencias para la evaluación: (\strong{holdout}) reticencia de datos son particionado en
#'  un conjunto de datos de entrenamiento y uno de evaluación usando una fracción (\emph{cv.param[1]})
#'  y la operación puede ser repetida \emph{cv.param[2]}) veces, \strong{LOO}
#'  (Leave One Out) cada punto es tomado sucesivamente como dato de evaluación.}
#'  \item{metric}{Choice of the metric used to compute the binary map threshold
#'  and the confusion matrix (by default SES as recommended by Liu et al.
#'  (2005), see reference below): \strong{Kappa} maximizes the Kappa,
#'  \strong{CCR} maximizes the proportion of correctly predicted observations,
#'  \strong{TSS} (True Skill Statistic) maximizes the sum of sensitivity and
#'  specificity, \strong{SES} uses the sensitivity-specificity equality,
#'  \strong{LW} uses the lowest occurrence prediction probability, \strong{ROC}
#'  minimizes the distance between the ROC plot (receiving operating
#'  characteristic curve) and the upper left corner (1,1).
#'  \strong{(Esp)} Elección de la métrica usada para procesar el umbral del mapa bianrio y la matriz de
#'  confusión (SES esta por defecto, recomentdado por Liu et al. (2005), ver referencias abajo):
#'  \strong{Kappa} maximiza el Kappa, \strong{CRR} maximiza la correcta proporción de las observaciones
#'  predecidas, \strong{TSS} (Verdaderas habilidades estaditicas) maximiza la suma de la especificidad
#'  y sensitividad, \strong{SES} usa la igualdad sensitividad-especificidad, \strong{LW} usa la
#'  menor probabilidad de predicción de ocurrencias, \strong{ROC} miniza la distancia entre la curva ROC
#'  (Características de funcionamiento del receptor) y la esquina superior izquierda  }\item{axes.metric}{Metric
#'  used to evaluate the variable relative importance (difference between a full model and one with
#'  each variable successively omitted): \strong{Pearson} (computes a simple Pearson's correlation \emph{r}
#'  between predictions of the full model and the one without a variable, and returns the score \emph{1-r}:
#'  the highest the value, the more influence the variable has on the model), \strong{AUC}, \strong{Kappa},
#'  \strong{sensitivity}, \strong{specificity}, and \strong{prop.correct} (proportion of correctly predicted
#'  occurrences).
#'  \strong{(Esp)} Métrica usada para evaluar la variable de importancia relativa (diferencia entre el modelo completo
#'  y un modelos con cada variable omitida): \strong{Pearson} (procesa una simple correlación \emph{r} entre las
#'  predicciones del modelo completo y el que no tiene una variable, regresa el puntaje \emph{1-r}:el mas alto
#'  valor representa a la variable con mas influencia en el modelo),\strong{AUC}, \strong{Kappa},
#'  \strong{sensitivity}, \strong{specificity}, y \strong{prop.correct} (proporción del correcto predictor de ocurrencias)}
#'  \item{ensemble.metric}{Ensemble metric(s) used to select SDMs: \strong{AUC}, \strong{Kappa}, \strong{sensitivity},
#'  \strong{specificity}, and \strong{prop.correct} (proportion of correctly predicted occurrences).}
#'  \strong{(Esp)} Usado para seleccionar del SDM el: \strong{AUC}, \strong{Kappa}, \strong{sensitivity},
#'  \strong{specificity}, y \strong{prop.correct} (proporción del correcto predictor de ocurrencias)
#'  \item{"..."}{See algorithm in detail section.
#'  \strong{(Esp)} Ver los algortimos en la sección de detalles}}
#'
#'@section Generalized linear model (\strong{GLM}) : Uses the \code{glm}
#'  function from the package 'stats', you can set the following parameters (see
#'  \code{\link[stats]{glm}} for more details): \describe{
#'  \item{test}{character. Test used to evaluate the SDM, default 'AIC'.}
#'  \item{epsilon}{numeric. Positive convergence tolerance eps ; the iterations
#'  converge when \emph{|dev - dev_{old}|/(|dev| + 0.1) < eps}. By default, set
#'  to 10e-08.} \item{maxit}{numeric. Integer giving the maximal number of IWLS
#'  (Iterative Weighted Last Squares) iterations, default 500.} }
#'
#'@section \strong{(Esp)} Modelo linear generalizado (\strong{GLM}) : Utiliza la función \code{glm} del paquete 'stats',
#'  se puede colocar los siguientes parametros (ver \code{\link[stats]{glm}} para mas detalles):
#'  \describe{\item{test}{cadena. Test usado para evaluar el SDM, por defecto 'AIC'.} \item{epsilon}{ número. Tolerancia
#'  de convergencia positiva eps ; las iteraciones convergen cuando \emph{|dev - dev_{old}|/(|dev| + 0.1) < eps}.
#'  por defecto, se establece en 10e-08.} \item{maxit}{número. Entero que representa el máximo número de IWLS
#'  (Últimos iterativos cuadrados ponderados) iteraciones, por defecto 500.} }
#'
#'@section Generalized additive model for large dataset(\strong{BAM}) : Uses the \code{bam}
#'  function from the package 'mgcv', you can set the following parameters (see
#'  \code{\link[mgcv]{bam}} for more details): \describe{ \item{test}{character.
#'  Test used to evaluate the model, default 'AIC'.} \item{epsilon}{numeric.
#'  This is used for judging conversion of the GLM IRLS (Iteratively Reweighted
#'  Least Squares) loop, default 10e-08.} \item{maxit}{numeric. Maximum number
#'  of IRLS iterations to perform, default 500.} }
#'
#'@section Modelos aditivos generalizados para grandes datos(\strong{BAM}) : Utiliza la función \code{bam} del paquete
#'  'mgcv', se puede colocar los siguientes parametros (ver \code{\link[mgcv]{bam}} para más detalles): \describe{
#'  \item{test}{cadena. Test usado para evaluar el modelo, por defecto 'AIC'. } \item{epsilon}{número.
#'  Esto se utiliza para juzgar la conversión del bucle GLM IRLS (Mínimos cuadrados reponderados reiterativamente),
#'  por defecto 10e-08.} \item{maxit}{número. Máximo número de iteraciones IRLS a realizar, por defecto 500.} }
#'
#'@section Multivariate adaptive regression splines (\strong{MARS}) : Uses the
#'  \code{earth} function from the package 'earth', you can set the following
#'  parameters (see \code{\link[earth]{earth}} for more details): \describe{
#'  \item{degree}{integer. Maximum degree of interaction (Friedman's mi) ; 1
#'  meaning build an additive model (i.e., no interaction terms). By default, set to 2.} }
#'
#'@section \strong{(Esp)} Reglas multivariadas de regresión adaptativa (\strong{MARS}) : Utiliza la función \code{earth} del paquete
#'  'earth', se puede colocar los siguientes parámetros (ver \code{\link[earth]{earth}} para mas detalles): \describe{
#'  \item{degree}{entero. Máximo grado de interacción (mi de Friedman); 1 significa construir un modelo aditivo
#'  (es decir, sin términos de interacción). De forma predeterminada, se establece en 2.} }
#'
#'@section Boosted regressions model (\strong{GBM}) : Uses the
#'  \code{GBM} function from the package 'gbm' you can set the following
#'  parameters (see \code{\link[gbm]{gbm}} for more details): \describe{
#'  \item{iterations}{integer. The total number of trees to fit. This is equivalent
#'  to the number of iterations and the number of basis functions in the
#'  additive expansion. By default, set to 250.}}
#'
#'@section \strong{(Esp)} Modelo aumentado de regresión generalizado (\strong{GBM} : Utiliza la función \code{GBM}) del paquete
#'  'gbm' se puede colocar los siguientes parámetros (ver \code{\link[gbm]{gbm}} para más detalles) : \describe{
#'  \item{iterations}{número. El número total de árboles para encajar. Esto es equivalente al número de iteraciones
#'  y al número de funciones básicas en la expansión aditiva. De forma predeterminada, se establece en 250.}}
#'
#'@section Classification tree analysis (\strong{CTA}) : Uses the \code{rpart}
#'  function from the package 'rpart', you can set the following parameters (see
#'  \code{\link[rpart]{rpart}} for more details): \describe{
#'  \item{final.leave}{integer. The minimum number of observations in any
#'  terminal node, default 1.} \item{algocv}{integer. Number of
#'  cross-validations, default 3.} }
#'
#'@section \strong{(Esp)} Análisis de clasificación de árboles (\strong{CTA}) : Utiliza la función \code{rpart} del paquete
#'  'rpart', se puede colocar los siguientes parámetros (ver \code{\link[rpart]{rpart}} para más detalles):
#'  \describe{ \item{final.leave}{número. El mínimo número de observaciones en cualquier nodo final, por defecto
#'  en 1.} \item{algocv}{número. Número de validaciones cruzadas, por defecto 3.} }
#'
#'@section Random Forest (\strong{RF}) : Uses the \code{randomForest} function
#'  from the package 'randomForest', you can set the following parameters (see
#'  \code{\link[randomForest]{randomForest}} for more details): \describe{
#'  \item{trees}{integer. Number of trees to grow. This should not be set to a
#'  too small number, to ensure that every input row gets predicted at least a
#'  few times. By default, set to 500.} \item{final.leave}{integer. Minimum
#'  size of terminal nodes. Setting this number larger causes smaller trees to
#'  be grown (and thus take less time). By default, set to 1.} }
#'
#'@section \strong{(Esp)} Árboles Aleatorios (\strong{RF}) : Utiliza la función \code{randomForest} del paquete
#'  'randomForest', se puede colocar los siguientes parámetros (ver \code{\link[randomForest]{randomForest}}
#'  para más detalles): \describe{ \item{trees}{número, Número de árboles a crecer. Este no debe ser muy pequeño
#'  , para garantizar que cada fila de entrada se predice al menos unas cuantas veces. De forma predeterminada,
#'  se establece en 250. }}
#'
#'@section Maximum Entropy (\strong{MAXNET}) : Uses the \code{maxnet} function
#'  from the package 'maxnet'. Maxent species distribution modeling using glmnet for model fitting (see
#'  \code{\link[maxnet]{maxnet}} for more details).
#'
#'@section \strong{(Esp)} Máxima entropía (\strong{MAXNET}) : Utiliza la función \code{maxnet} del paquete 'maxnet'. Maxent modelamiento
#'  de distribución de especies que utiliza glmnet para el ajuste del modelo (ver \code{\link[maxnet]{maxnet}} para mas detalles)
#'
#'@section Artificial Neural Network (\strong{ANN}) : Uses the \code{nnet}
#'  function from the package 'nnet', you can set the following parameters (see
#'  \code{\link[nnet]{nnet}} for more details): \describe{ \item{maxit}{integer.
#'  Maximum number of iterations, default 500.} }
#'
#'@section \strong{(Esp)} Red neuronal artificial (\strong{ANN}) : Utiliza la función \code{nnet} del paquete 'nnet', se puede colocar los
#'siguientes parámetro (ver \code{\link[nnet]{nnet}} para más detalles): \describe{ \item{maxit}{
#'número. Número máximo de iteraciones, por defecto 500.} }
#'
#'@section Support vector machines (\strong{KSVM}) : Uses the \code{ksvm} function
#'  from the package 'kernlab', you can set the following parameters (see
#'  \code{\link[kernlab]{ksvm}} for more details): \describe{ \item{epsilon}{float.
#'  Epsilon parameter in the insensitive loss function, default 1e-08.}
#'  \item{algocv}{integer. If an integer value k>0 is specified, a k-fold
#'  cross-validation on the training data is performed to assess the quality of
#'  the model: the accuracy rate for classification and the Mean Squared Error
#'  for regression. By default, set to 3.} }
#'
#'@section \strong{(Esp)} Máquinas vectoriales de apoyo (\strong{KSVM}) : Utiliza la función del paquete 'kernlab', se
#'  puede colocar los siguientes parámetros (ver \code{\link[kernlab]{ksvm}} para más detalles) : \describe{
#'  \item{epsilon}{Número flotante. Parámetro de Epsilon en la función de pérdida insensible, por defecto 1e-08.}
#'  \item{algocv}{número. Si se especifica un valor entero k> 0, se realiza una validación cruzada de k veces en
#'  los datos de entrenamiento para evaluar la calidad del modelo: la tasa de precisión para la clasificación y
#'  el error cuadrático medio para la regresión. De forma predeterminada, se establece en 3.}}
#'
#'@section Warning : Depending on the raster object resolution the process can
#'  be more or less time and memory consuming.
#'\strong{(Esp)} Advertencia : dependiendo da la resolución del objeto ráster el proceso
#' puede tardar mas o menos tiempo y consumo de memoria.
#'
#' @examples
#' # Loading data
#' data(Env)
#' data(Occurrences)
#' Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
#'
#' # SDM building
#' SDM <- modelling('GLM', Occurrences, Env, Xcol = 'LONGITUDE', Ycol = 'LATITUDE')
#'
#' # Results plotting
#' \dontrun{
#' plot(SDM)
#' }
#'
#'
#'@seealso \code{\link{ensemble_modelling}} to build ensemble SDMs,
#'  \code{\link{stack_modelling}} to build SSDMs.
#'
#'@references M. Barbet-Massin, F. Jiguet, C. H.  Albert, & W. Thuiller (2012)
#'  'Selecting pseudo-absences for species distribution models: how, where and
#'  how many?' \emph{Methods Ecology and Evolution} 3:327-338
#'  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x/full}
#'
#'
#'
#'
#'
#'
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) 'Selecting
#'  thresholds of occurrence in the prediction of species distributions.'
#'  \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'
#'
#'
#'
#'
#'
#'@export
modelling <- function(algorithm, Occurrences, Env, Xcol = "Longitude",
                      Ycol = "Latitude", Pcol = NULL, name = NULL, PA = NULL, cv = "holdout",
                      cv.param = c(0.7, 2), thresh = 1001, metric = "SES", axes.metric = "Pearson",
                      select = FALSE, select.metric = c("AUC"), select.thresh = c(0.75),
                      verbose = TRUE, GUI = FALSE, ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, name = name, PA = PA,
             cv = cv, cv.param = cv.param, thresh = thresh, metric = metric,
             axes.metric = axes.metric, select = select, select.metric = select.metric,
             select.thresh = select.thresh, verbose = verbose, GUI = GUI)

  # Test if algorithm is available
  available.algo <- c("GLM", "BAM", "MARS", "GBM", "CTA", "RF", "MAXNET",
                      "ANN", "KSVM")
  if (algorithm == "all") {
    algorithm <- available.algo
  }
  if (!(algorithm %in% available.algo)) {
    stop(algorithm, " is still not available, please use one of those : GLM, BAM, MARS, GBM, CTA, RF, MAXNET, ANN, KSVM")
  }

  # Empty Algorithm niche model object creation
  model <- Algorithm.SDM(algorithm)
  if (!is.null(name)) {
    name <- paste0(name, ".")
  }
  model@name <- paste0(name, algorithm, ".SDM")
  model@parameters$data <- "presence/absence data set"
  model@parameters$metric <- metric

  if (verbose) {
    cat("Data check ... \n")
  }
  # Occurrences data input test | Data frame needed
  if (is.matrix(Occurrences)) {
    Occurrences <- data.frame(Occurrences)
  }
  if (!is.data.frame(Occurrences)) {
    stop("Occurrences data set is not a data frame or a matrix")
  }
  if ((Xcol %in% names(Occurrences)) == FALSE) {
    stop("X column is not well defined")
  }
  if ((Ycol %in% names(Occurrences)) == FALSE) {
    stop("Y column is not well defined")
  }
  if (is.null(Pcol)) {
    PO <- TRUE  # Presence only
    if (verbose) {
      cat("No presence column, presence-only data set is supposed.\n")
    }
    model@parameters$data <- "presence-only data set"
  } else if ((Pcol %in% names(Occurrences)) == FALSE) {
    stop("Presence column is not well defined")
  } else {
    PO <- FALSE
  }
  if (!is.null(PA)) {
    PO <- TRUE
  }
  if (PO) {
    if (verbose) {
      cat("Pseudo-absence selection will be computed.\n")
    }
  }
  data <- data.frame(X = Occurrences[which(names(Occurrences) == Xcol)],
                     Y = Occurrences[which(names(Occurrences) == Ycol)])
  names(data) <- c("X", "Y")
  if (PO) {
    data$Presence <- 1
  } else {
    data$Presence <- Occurrences[, which(names(Occurrences) == Pcol)]
  }

  # Environment data input test | RasterStack needed
  if (inherits(Env, "Raster")) {
    Env <- stack(Env)
  }
  if (!inherits(Env, "RasterStack")) {
    stop("Environment data set is not a raster or a raster stack")
  }
  if (verbose) {
    cat("   done. \n\n")
  }
  if (GUI) {
    incProgress(1/5, detail = "Data ckecked")
  }

  # Pseudo - absences selection
  model@data <- data
  if (PO) {
    if (verbose) {
      cat("Pseudo absence selection... \n")
    }
    model <- PA.select(model, Env, PA, verbose)
    model@parameters["PA"] <- TRUE
    if (verbose) {
      cat("   done. \n\n")
    }
    if (GUI) {
      incProgress(1/5, detail = "Pseudo-absence selected")
    }
  }
  model <- data.values(model, Env)

  # Evaluation
  if (verbose) {
    cat("Model evaluation...\n")
  }
  model <- evaluate(model, cv, cv.param, thresh, metric, Env, ...)
  if (verbose) {
    cat("   done. \n\n")
  }
  if (GUI) {
    incProgress(1/5, detail = "SDM evaluated")
  }

  # Model selection
  test <- TRUE
  if (select) {
    for (j in seq_len(length(select.metric))) {
      if (model@evaluation[, which(names(model@evaluation) == select.metric[j])] <
          select.thresh[j]) {
        test <- FALSE
      }
    }
  }
  if (test) {
    # Projection
    if (verbose) {
      cat("Model projection...\n")
    }
    model <- project(model, Env, ...)
    if (verbose) {
      cat("   done. \n\n")
    }
    if (GUI) {
      incProgress(1/5, detail = "SDM projected")
    }

    # Axes evaluation
    if (verbose) {
      cat("Model axes contribution evaluation...\n")
    }
    model <- evaluate.axes(model, cv, cv.param, thresh, metric, axes.metric,
                           Env, ...)
    if (verbose) {
      cat("   done. \n\n")
    }
    if (GUI) {
      incProgress(1/5, detail = "SDM axes contribution evaluated")
    }
    rm(list = ls()[-which(ls() == "model")])
    gc()
    return(model)
  } else {
    if (verbose) {
      cat("Model have been rejected, NULL is returned ! \n")
    }
    if (GUI) {
      incProgress(2/5, detail = "SDM rejected")
    }
    rm(list = ls())
    gc()
    return(NULL)
  }
}
