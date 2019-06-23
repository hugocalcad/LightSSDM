#' @include Algorithm.SDM.R ensemble.R Ensemble.SDM.R checkargs.R
#'   stacking.R Stacked.SDM.R
#' @importFrom shiny incProgress
#' @importFrom raster stack writeRaster
NULL

#'Build an SSDM that assembles multiple algorithms and species.
#'\strong{(Esp)} Construye un SSDM que ensambla mutilples algoritmos y especies.
#'
#'This is a function to build an SSDM that assembles multiple algorithm and
#'species. The function takes as inputs an occurrence data frame made of
#'presence/absence or presence-only records and a raster object for data
#'extraction and projection. The function returns an S4
#'\linkS4class{Stacked.SDM} class object containing the local species richness
#'map, the between-algorithm variance map, and all evaluation tables coming with
#'(model evaluation, algorithm evaluation, algorithm correlation matrix and
#'variable importance), and a list of ensemble SDMs for each species (see
#'\code{\link{ensemble_modelling}}).
#'\strong{(Esp)} Esta función para construir un SSDM que ensambla mutiples algoritmos y especies.
#'La toma como entradas las ocurrencias de presencias/ausencias o precencia solamente y
#'un objeto ráster para extracción y proyección de datos. La función retorna objeto de clase
#'S4 \linkS4class{Stacked.SDM} que contiene el mapa de riqueza de especies, el mapa de
#'varianza entre algoritmos, y todos las tablas de evaluación que contiene (evaluación del
#'modelo, evaluación de algoritmo, matriz de correlación de algoritmos y variable de
#'importancia), y un lista de SDMs ensamblados para cada especie (ver
#'\code{\link{ensemble_modelling}}).
#'
#'@param algorithms character. Choice of the algorithm(s) to be run (see details below).
#'\strong{(Esp)} Elección de los algoritmos que se ejecutarán (ver detalles más abajo).
#'@param Occurrences data frame. Occurrence table (can be processed first by \code{\link{load_occ}}).
#'\strong{(Esp)} Tabla de ocurrencias (puede ser procesado primero por \code{\link{load_occ}}).
#'@param Env raster object. Raster object of environmental variables (can be
#'  processed first by \code{\link{load_var}}).
#'\strong{(Esp)} Objeto ráster de las variables ambientales (puede ser procesado primero
#'  por \code{\link{load_var}}).
#'@param Xcol character. Name of the column in the occurrence table containing Latitude or X coordinates.
#'\strong{(Esp)} Nombre de la columna en la tabla de ocurrencias que contiene la latitud o la coordenada X.
#'@param Ycol character. Name of the column in the occurrence table containing Longitude or Y coordinates.
#'\strong{(Esp)} Nombre de la columna en la tabla de ocurrencias que contiene la longitud o la coordenada Y.
#'@param Pcol character. Name of the column in the occurrence table specifying
#'  whether a line is a presence or an absence. A value of 1 is presence and
#'  value of 0 is absence. If NULL presence-only dataset is assumed.
#'\strong{(Esp)} Nombre de la columna en la tabla de ocurrencias especificando
#'  si la linea es una presencia o ausencia. Un valor de 1 es presencia y un valor
#'  de 0 es ausencia. Si es NULL se asume solo datos de presencia.
#'@param Spcol character. Name of the column containing species names or IDs.
#'\strong{(Esp)} Nombre de la columna que contien el nombre de la especie o los IDs.
#'@param rep integer. Number of repetitions for each algorithm.
#'\strong{(Esp)} Número de repeticiones para cada algortimo.
#'@param name character. Optional name given to the final Ensemble.SDM produced.
#'\strong{(Esp)} Nombre opcional dado al Ensemble.SDM final producido.
#'@param save logical. If set to true, the SSDM is automatically saved.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}), el SSDM es automaticamente guardado.
#'@param path character. If save is true, the path to the directory in which the ensemble SDM will be saved.
#'\strong{(Esp)} Si guardar (save) es verdadero (\code{TRUE}), la ruta al directorio donde el SDM ensamblado sera guardado.
#'@param PA list(nb, strat) defining the pseudo-absence selection strategy used
#'  in case of presence-only dataset. If PA is NULL, recommended PA selection
#'  strategy is used depending on the algorithm (see details below).
#'\strong{(Esp)} definieniendo la estrategia de selección de las pseudo-ausencias usadas
#'  en caso de datos con solo presencias. Si PA es nulo (NULL), estrategias recomendadas
#'  son usadas dependiendo del algoritmo (ver detalles más abajo).
#'@param cv character. Method of cross-validation used to evaluate the ensemble SDM (see details below).
#'\strong{(Esp)} Metodo de validación-cruzada usada para evaluar el SDM ensamblado (ver detalles más abajo).
#'@param cv.param numeric. Parameters associated to the method of
#'  cross-validation used to evaluate the ensemble SDM (see details below).
#'\strong{(Esp)} Parametros asociados al método de validación- cruzada usada
#'  para evaluar el SDM ensamblado (ver detalles más abajo).
#'@param thresh numeric. A single integer value representing the number of equal
#'  interval threshold values between 0 and 1 (see \code{\link[SDMTools]{optim.thresh}}).
#'\strong{(Esp)} Un valor numérico simple  representado el numero de igual intervalo
#'  del umbral entre 0 y 1 (ver \code{\link[SDMTools]{optim.thresh}}).
#'@param metric character. Metric used to compute the binary map threshold (see details below.)
#'\strong{(Esp)} Métrica usada para procesar el umbral del mapa binario(ver detalles más abajo)
#'@param axes.metric Metric used to evaluate variable relative importance (see details below).
#'\strong{(Esp)} Métrica usada para evaluar la variable relativa de importancia (ver detalles más abajo)
#'@param uncertainty logical. If \code{TRUE}, generates an uncertainty map and an algorithm correlation matrix.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}), se genera un mapa de incertidumbre y una matriz de correlación de algoritmos.
#'@param tmp logical. If set to true, the habitat suitability map of each
#'  algorithm is saved in a temporary file to release memory. But beware: if you
#'  close R, temporary files will be deleted To avoid any loss you can save your
#'  ensemble SDM with \code{\link{save.model}}. Depending on number, resolution
#'  and extent of models, temporary files can take a lot of disk space.
#'  Temporary files are written in R environment temporary folder.
#'\strong{(Esp)} Si es verdadero (true), el mapa de habitabilidad de cada algortimo es guardado
#'  en una carpeta temporal para liberar memoria. Pero hay que tener cuidado: si cierras
#'  R,los archivos temporales se´ran eliminadospara evitar cualquier perdida  puedes guardar
#'  tus SDMs ensamblados con \code{\link{save.model}}.  dependiendo al numero, resolución
#'  y extención de modelos,los archivos temporales pueden ocupar mucho espacio en disco.
#'  Los archivos temporales son escritos en el contexto de archivos temporales de R.
#'@param ensemble.metric character. Metric(s) used to select the best SDMs that
#'  will be included in the ensemble SDM (see details below).
#'\strong{(Esp)} Métrica(s) usada para seleccionar el mejor SDMs que puede ser incluido en el
#'  SDM ensamblado (ver detalles más abajo).
#'@param ensemble.thresh numeric. Threshold(s) associated with the metric(s)
#'  used to compute the selection.
#'\strong{(Esp)} Umbral(es) asociado con la métrica(s) usado para procesar la selección
#'@param weight logical. If \code{TRUE}, SDMs are weighted using the ensemble
#'  metric or, alternatively, the mean of the selection metrics.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}), SDMs son calificadas usando la métrica de
#'  ensamblado o, alternativamente, el promedio de las métricas de selección.
#'@param method character. Define the method used to create the local species
#'  richness map (see details below).
#'\strong{(Esp)} define el método usado para crear el mapa local de riqueza de especies
#'  (ver detalles más abajo).
#'@param rep.B integer. If the method used to create the local species richness
#'  is the random bernoulli (\strong{Bernoulli}), rep.B parameter defines the number of
#'  repetitions used to create binary maps for each species.
#'\strong{(Esp)} Si el método usado para crear el mapa local de riqueza de especies es el
#'  Bernoulli aleatorio (\strong{Bernoulli}), el parametro rep.B  define el número de
#'  repeticiones usados para crear el mapa binario para cada especie.
#'@param range integer. Set a value of range restriction (in pixels) around
#'  presences occurrences on habitat suitability maps (all further points will
#'  have a null probability, see Crisp et al (2011) in references). If NULL, no
#'  range restriction will be applied.
#'\strong{(Esp)} Establece un valor de un rango de restricción (en pixeles) cerca a la
#' presencia de ocurrencias en el mapa de idoneidad de hábitat (todos los puntos lejanos
#' tendran una probabilidad nula, ver Crisp et al (2011) en referencias). Si es NULL, no
#' se aplicara ninguna restricción de rango.
#'@param endemism character. Define the method used to create an endemism map (see details below).
#'\strong{(Esp)} define el metodo usado para crear un mapa de endemismo (ver detalles más abajo).
#'@param verbose logical. If set to true, allows the function to print text in the console.
#'\strong{(Esp)} Si es verdadero (\code{TRUE}), permite la función de imprimir texto en la consola.
#'@param GUI logical. Don't take that argument into account (parameter for the user interface).
#'\strong{(Esp)} No tomar este arguento en cuenta(Parametro para la interfaz de usuario)
#'@param cores integer. Specify the number of CPU cores used to do the
#'  computing. You can use \code{\link[parallel]{detectCores}}) to automatically
#'  used all the available CPU cores.
#'\strong{(Esp)} Especifica el número de cores del CPU usados para el procesamiento. puedes usar
#'  \code{\link[parallel]{detectCores}})  para usar automaticamente todos los cores disponibles del CPU.
#'@param ... additional parameters for the algorithm modelling function (see details below).
#'\strong{(Esp)} parametros adicionales la funcion modelling de algoritmos (ver detalles más abajo)
#'
#'@return an S4 \linkS4class{Stacked.SDM} class object viewable with the \code{\link{plot.model}} function.
#'\strong{(Esp)} un objeto de clase \linkS4class{Stacked.SDM} visible con la función \code{\link{plot.model}}
#'
#'@details \describe{ \item{algorithms}{'all' calls all the following
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
#'  \strong{(Esp)} Lista con dos valores: \strong{nb} número de pseudo-ausencias seleccionados, y \strong{strat}
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
#'  (Características de funcionamiento del receptor) y la esquina superior izquierda. }
#'  \item{axes.metric}{Metric used to evaluate the variable relative importance (difference between a full model and one with
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
#'  \strong{(Esp)} Usado para seleccionar del SDM el: \strong{AUC}, \strong{Kappa}, \strong{sensitivity,
#'  \strong{specificity}, y \strong{prop.correct} (proporción del correcto predictor de ocurrencias)}
#'  \item{method}{Choice of the method used to compute the local
#'  species richness map (see Calabrese et al. (2014) and D'Amen et al (2015) for
#'  more informations, see reference below): \strong{(Esp)} Elige el método usado para procesar
#'  el mapa de riqueza de especies (ver Calabrese et al. (2014) y D'Amen el al (2015) para
#'  más información ver referencias abajo): \strong{pSSDM} sum
#'  probabilities of habitat suitability maps. \strong{(Esp)} suma de probabilidades de los
#'  mapas de idoneidad del hábitat. \strong{Bernoulli} draw repeatedly
#'  from a Bernoulli distribution. \strong{(Esp)} dibujar repetidamente una distribución de Bernoulli.
#'  \strong{bSSDM} sum the binary map obtained with the thresholding (depending on the metric of the
#'  ESDM). \strong{(Esp)} suma el mapa binario obtenido con el umbral (según la métrica del ESDM).
#'  \strong{MaximumLikelihood} adjust species richness of the model by
#'  linear regression. \strong{(Esp)} Ajustar la riqueza de especies del modelo mediante regresión lineal.
#'  \strong{PRR.MEM} model richness with a macroecological model (MEM) and adjust each ESDM binary map by
#'  ranking habitat suitability and keeping as much as predicted richness of the MEM. \strong{(Esp)}
#'  Modela la riqueza con un modelo macroecológico (MEM) y ajuste cada mapa binario de ESDM clasificando
#'  la aptitud del hábitat y manteniendo la riqueza predicha del MEM. \strong{PRR.pSSDM} model
#'  richness with a pSSDM and adjust each ESDM binary map by ranking habitat
#'  suitability and keeping as much as predicted richness of the pSSDM. \strong{(Esp)} Modela la riqueza con un
#'  pSSDM y ajuste cada mapa binario de ESDM clasificando la idoneidad del hábitat y manteniendo la riqueza predicha del pSSDM}
#'  \item{endemism}{Choice of the method used to compute the endemism map (see
#'  Crisp et al. (2001) for more information, see reference below):
#'  \strong{NULL} No endemism map, \strong{WEI} (Weighted Endemism Index)
#'  Endemism map built by counting all species in each cell and weighting each
#'  by the inverse of its range, \strong{CWEI} (Corrected Weighted Endemism
#'  Index) Endemism map built by dividing the weighted endemism index by the
#'  total count of species in the cell. First string of the character is the
#'  method either WEI or CWEI, and in those cases second string of the vector is
#'  used to precise range calculation, whether the total number of occurrences
#'  \strong{'NbOcc'} whether the surface of the binary map species distribution
#'  \strong{'Binary'}.
#'  \strong{(Esp)} Elige el método usado para procesar el mapa de endemismo  (ver
#'  Crisp et al. (2001) para mas infformación, ver referencias más abajo):
#'  \strong{NULL} Sin mapa de endemismo, \strong{WEI} (Indice de endemismo ponderado)
#'  Mapa de endemismo construido al contar todas las especies en cada celda y ponderar
#'  cada una por el inverso de su rango, \strong{CWEI} (Indice de endemismo ponderado corregido)
#'  el mapa de endemismo construido dividiendo el índice de endemismo ponderado por el recuento total
#'  de especies en la celda. La primera cadena del carácter es el método WEI o CWEI, y en esos casos
#'  la segunda cadena del vector se usa para el cálculo preciso del rango, ya sea el número total de ocurrencias
#'  \strong{'NbOcc'} si la superficie de la distribución de las especies del mapa binario \strong{'Binary'}.}
#'  \item{...}{See algorithm in detail section. \strong{(Esp)} ver lagortimo en la sección de detalles.} }
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
#'  se puede colocar los siguientes parametros (ver \code{\link[stats]{glm}} para mas detalles): \describe{
#'  \item{test}{cadena. Test usado para evaluar el SDM, por defecto 'AIC'.} \item{epsilon}{ número. Tolerancia
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
#'  parameters (see \code{\link[earth]{earth}} for more details): \describe{\item{degree}{integer.
#'  Maximum degree of interaction (Friedman's mi) ; 1
#'  meaning build an additive model (i.e., no interaction terms). By default, set to 2.} }
#'
#'@section \strong{(Esp)} Reglas multivariadas de regresión adaptativa (\strong{MARS}) : Utiliza la función \code{earth} del paquete
#'  'earth', se puede colocar los siguientes parámetros (ver \code{\link[earth]{earth}} para mas detalles):
#'  \describe{\item{degree}{entero. Máximo grado de interacción (mi de Friedman); 1 significa construir un modelo aditivo
#'  (es decir, sin términos de interacción). De forma predeterminada, se establece en 2.} }
#'
#'@section Boosted regressions model (\strong{GBM}) : Uses the
#'  \code{GBM} function from the package 'gbm' you can set the following
#'  parameters (see \code{\link[gbm]{gbm}} for more details): \describe{\item{iterations}{integer.
#'  The total number of trees to fit. This is equivalent
#'  to the number of iterations and the number of basis functions in the
#'  additive expansion. By default, set to 250.}}
#'
#'@section \strong{(Esp)} Modelo aumentado de regresión generalizado (\strong{GBM} : Utiliza la función \code{GBM}) del paquete
#'  'gbm' se puede colocar los siguientes parámetros (ver \code{\link[gbm]{gbm}} para más detalles) :
#'  \describe{\item{iterations}{número. El número total de árboles para encajar. Esto es equivalente al número de iteraciones
#'  y al número de funciones básicas en la expansión aditiva. De forma predeterminada, se establece en 250.}}
#'
#'@section Classification tree analysis (\strong{CTA}) : Uses the \code{rpart}
#'  function from the package 'rpart', you can set the following parameters (see
#'  \code{\link[rpart]{rpart}} for more details): \describe{\item{final.leave}{integer.
#'  The minimum number of observations in any
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
#'  \code{\link[randomForest]{randomForest}} for more details):
#'  \describe{\item{trees}{integer. Number of trees to grow. This should not be set to a
#'  too small number, to ensure that every input row gets predicted at least a
#'  few times. By default, set to 500.} \item{final.leave}{integer. Minimum
#'  size of terminal nodes. Setting this number larger causes smaller trees to
#'  be grown (and thus take less time). By default, set to 1.} }
#'
#'@section \strong{(Esp)} Árboles Aleatorios (\strong{RF}) : Utiliza la función \code{randomForest} del paquete
#'  'randomForest', se puede colocar los siguientes parámetros (ver \code{\link[randomForest]{randomForest}}
#'  para más detalles): \describe{\item{trees}{número, Número de árboles a crecer. Este no debe ser muy pequeño
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
#'  \code{\link[nnet]{nnet}} for more details): \describe{\item{maxit}{integer.
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
#'
#'@section \strong{(Esp)} Advertencia : dependiendo da la resolución del objeto ráster el proceso
#' puede tardar mas o menos tiempo y consumo de memoria.
#'
#' @examples
#' \dontrun{
#' # Loading data
#' data(Env)
#' data(Occurrences)
#'
#' # SSDM building
#' SSDM <- stack_modelling(c('CTA', 'KSVM'), Occurrences, Env, rep = 1,
#'                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
#'                        Spcol = 'SPECIES')
#'
#' # Results plotting
#' plot(SSDM)
#' }
#'
#'@seealso \code{\link{modelling}} to build simple SDMs.
#'
#'
#'@references M. D'Amen, A. Dubuis, R. F. Fernandes, J. Pottier, L. Pelissier, &
#'  A Guisan (2015) "Using species richness and functional traits prediction to
#'  constrain assemblage predicitions from stacked species distribution models"
#'  \emph{Journal of Biogeography} 42(7):1255-1266
#'  \url{http://doc.rero.ch/record/235561/files/pel_usr.pdf}
#'
#'  M. Barbet-Massin, F. Jiguet, C. H.  Albert, & W. Thuiller (2012) "Selecting
#'  pseudo-absences for species distribution models: how, where and how many?"
#'  \emph{Methods Ecology and Evolution} 3:327-338
#'  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x/full}
#'
#'
#'
#'  J.M. Calabrese, G. Certain, C.  Kraan, & C.F. Dormann (2014) "Stacking
#'  species distribution  models  and  adjusting  bias  by linking them to
#'  macroecological models." \emph{Global Ecology and Biogeography} 23:99-112
#'  \url{http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/calabrese2013globalecolbiogeogr.pdf}
#'
#'
#'
#'  M. D. Crisp, S. Laffan, H. P. Linder & A. Monro (2001) "Endemism in the
#'  Australian flora"  \emph{Journal of Biogeography} 28:183-198
#'  \url{http://biology-assets.anu.edu.au/hosted_sites/Crisp/pdfs/Crisp2001_endemism.pdf}
#'
#'
#'
#'  C. Liu, P. M. Berry, T. P. Dawson,  R. & G. Pearson (2005) "Selecting
#'  thresholds of occurrence in the prediction of species distributions."
#'  \emph{Ecography} 28:85-393
#'  \url{http://www.researchgate.net/publication/230246974_Selecting_Thresholds_of_Occurrence_in_the_Prediction_of_Species_Distributions}
#'
#'
#'@export
stack_modelling <- function(algorithms,
                           # Modelling data input
                           Occurrences, Env,
                           # Occurrences reading
                           Xcol = 'Longitude', Ycol = 'Latitude', Pcol = NULL, Spcol = 'SpeciesID',
                           # Model creation
                           rep = 10, name = NULL, save = FALSE, path = getwd(),
                           # Pseudo-absences definition
                           PA = NULL,
                           # Evaluation parameters
                           cv = 'holdout', cv.param = c(0.7,1), thresh = 1001,
                           axes.metric = 'Pearson', uncertainty = TRUE, tmp = FALSE,
                           # Assembling parameters
                           ensemble.metric = c('AUC'), ensemble.thresh = c(0.75), weight = TRUE,
                           # Diversity map computing
                           method = 'pSSDM', metric = 'SES', rep.B = 1000,
                           # Range restriction and endemism
                           range = NULL, endemism = c('WEI','Binary'),
                           # Informations parameters
                           verbose = TRUE, GUI = FALSE, cores = 1,
                           # Modelling parameters
                           ...) {
  # Check arguments
  .checkargs(Xcol = Xcol, Ycol = Ycol, Pcol = Pcol, Spcol = Spcol, rep = rep,
             name = name, save = save, path = path, PA = PA, cv = cv, cv.param = cv.param,
             thresh = thresh, axes.metric = axes.metric, uncertainty = uncertainty,
             tmp = tmp, ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
             weight = weight, method = method, metric = metric, rep.B = rep.B, range = range,
             endemism = endemism, verbose = verbose, GUI = GUI, cores = cores)

  # Test if algorithm is available
  available.algo <- c("GLM", "BAM", "MARS", "GBM", "CTA", "RF", "MAXNET",
                      "ANN", "KSVM")
  if ("all" %in% algorithms) {
    algorithms <- available.algo
  }
  for (i in seq_len(length(algorithms))) {
    if (!(algorithms[[i]] %in% available.algo)) {
      stop(algorithms[[i]], " is still not available, please use one of those : GLM, BAM, MARS, GBM, CTA, RF, MAXNET, ANN, KSVM")
    }
  }
  if (tmp) {
    cat("temporal..... \n")
    tmppath <- get("tmpdir", envir = .PkgEnv)
    if (!("/.enms" %in% list.dirs(tmppath)))
      (dir.create(paste0(tmppath, "/.enms")))
  }

  # Ensemble models creation
  if (verbose) {
    cat("#### Ensemble models creation ##### \n\n")
  }
  species <- levels(as.factor(Occurrences[, which(names(Occurrences) == Spcol)]))

  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (verbose) {
      cat("Opening clusters,", cores, "cores \n")
    }
    if ((parallel::detectCores() - 1) < cores) {
      warning("It seems you attributed more cores than your CPU have !")
    }
    cl <- parallel::makeCluster(cores, outfile = "")
    if (verbose) {
      cat("Exporting environment to clusters \n")
    }
    parallel::clusterExport(cl, varlist = c(lsf.str(envir = globalenv()),
                                            ls(envir = environment())), envir = environment())
    enms <- parallel::parLapply(cl, species, function(species) {
      enm.name <- species
      Spoccurrences <- subset(Occurrences, Occurrences[which(names(Occurrences) ==
                                                               Spcol)] == species)
      if (verbose) {
        cat("Ensemble modelling :", enm.name, "\n\n")
      }
      enm <- try(ensemble_modelling(algorithms, Spoccurrences, Env, Xcol,
                                    Ycol, Pcol, rep = rep, name = enm.name, save = FALSE, path = path,
                                    PA = PA, cv = cv, cv.param = cv.param, thresh = thresh, metric = metric,
                                    axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
                                    ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
                                    weight = weight, verbose = verbose, GUI = FALSE, n.cores = 1,
                                    ...))
      if (GUI) {
        incProgress(1/(length(levels(as.factor(Occurrences[, which(names(Occurrences) ==
                                                                     Spcol)]))) + 1), detail = paste(species, " ensemble SDM built"))
      }
      if (inherits(enm, "try-error")) {
        if (verbose) {
          cat(enm)
        }
        enm <- NULL
      } else {
        if (tmp && !is.null(enm)) {
          enm@projection <- writeRaster(enm@projection[[1]], paste0(tmppath,
                                                                    "/.enms/proba", enm.name), overwrite = TRUE)
          enm@binary <- writeRaster(enm@binary[[1]], paste0(tmppath,
                                                            "/.enms/bin", enm.name), overwrite = TRUE)
          enm@uncertainty <- writeRaster(enm@uncertainty, paste0(tmppath,
                                                                 "/.enms/uncert", enm.name), overwrite = TRUE)
        }
        if (verbose) {
          cat("\n\n")
        }
      }

      return(enm)
    })
    if (verbose) {
      cat("Closing clusters \n")
    }
    parallel::stopCluster(cl)

  } else {
    enms <- lapply(species, function(species) {
      enm.name <- species
      Spoccurrences <- subset(Occurrences, Occurrences[which(names(Occurrences) ==
                                                               Spcol)] == species)
      if (verbose) {
        cat("Ensemble modelling :", enm.name, "\n\n")
      }
      enm <- try(ensemble_modelling(algorithms, Spoccurrences, Env, Xcol,
                                    Ycol, Pcol, rep = rep, name = enm.name, save = FALSE, path = path,
                                    PA = PA, cv = cv, cv.param = cv.param, thresh = thresh, metric = metric,
                                    axes.metric = axes.metric, uncertainty = uncertainty, tmp = tmp,
                                    ensemble.metric = ensemble.metric, ensemble.thresh = ensemble.thresh,
                                    weight = weight, verbose = verbose, GUI = FALSE, ...))
      if (GUI) {
        incProgress(1/(length(levels(as.factor(Occurrences[, which(names(Occurrences) ==
                                                                     Spcol)]))) + 1), detail = paste(species, " ensemble SDM built"))
      }
      if (inherits(enm, "try-error")) {
        if (verbose) {
          cat(enm)
        }
        enm <- NULL
      } else {
        if (tmp && !is.null(enm)) {
          enm@projection <- writeRaster(enm@projection[[1]], paste0(tmppath,
                                                                    "/.enms/proba", enm.name), overwrite = TRUE)
          enm@uncertainty <- writeRaster(enm@uncertainty, paste0(tmppath,
                                                                 "/.enms/uncert", enm.name), overwrite = TRUE)
        }
        if (verbose) {
          cat("Terminado ",enm.name,"\n\n")
        }
      }
      return(enm)
    })
  }

  enms <- enms[!sapply(enms, is.null)]

  # Species stacking
  if (length(enms) < 2) {
    if (verbose) {
      stop("Less than two species models were retained, you should lower the ensemble threshold value (ensemble.thresh parameter).")
    } else {
      return(NULL)
    }
  } else {
    if (verbose) {
      cat("#### Species stacking with ensemble models ##### \n\n")
    }
    if (!is.null(name)) {
      enms["name"] <- name
    }
    enms["method"] <- method
    enms["rep.B"] <- rep.B
    if (method %in% c("PRR.MEM", "PRR.pSSDM")) {
      enms["Env"] <- Env
    }
    if (!is.null(range)) {
      enms["range"] <- range
    }
    enms$endemism <- endemism
    enms["verbose"] <- verbose
    stack <- do.call(stacking, enms)
  }

  if (!is.null(stack)) {
    # Paremeters
    stack@parameters$sp.nb.origin <- length(levels(as.factor(Occurrences[,
                                                                         which(names(Occurrences) == Spcol)])))
    if (GUI) {
      incProgress(1/(length(levels(as.factor(Occurrences[, which(names(Occurrences) ==
                                                                   Spcol)]))) + 1), detail = "SSDM built")
    }

    # Saving
    if (save) {
      if (verbose) {
        cat("#### Saving ##### \n\n")
      }
      if (!is.null(name)) {
        save.stack(stack, name = name, path = path)
      } else {
        save.stack(stack, path = path)
      }
    }
  }
  return(stack)
}
