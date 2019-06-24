#' LigthSSDM: Stacked species distribution modelling for large dataset
#' \strong{(Esp)} LightSSDM:  Apilado de modelos de distribución de especies (SSDM) para conjunto de datos grandes
#'
#' LightSSDM is a package to map species richness and endemism based on SSDM packages created by
#' Sylvain Schmitt (\url{https://github.com/sylvainschmitt/SSDM}). This package has the same
#' purpose as the original, but some algorithms were changed and some parameters were
#' changed to make it lighter and so we could work with more extensive data, instead of
#' the GAM function we used the BAM function of the mgcv package, instead of the MAXENT
#' function of the Dismo package we used the maxnet function of the maxnet package created
#' by the same developers of the MAXENT made in JAVA, instead of the svm function
#' of the e1071 package we used the ksvm function of the kernlab package.
#' \strong{(Esp)} LightSSDM es un paquete para mapear la riqueza y endemismo de especies basado en el paquete SSDM
#' creado por Sylvain Schmitt (\url{https://github.com/sylvainschmitt/SSDM}). Este paquete
#' tiene el mismo propósito que el original, peros se cambiaron algunos algoritmos y a otros
#' se le cambiaron algunos parámetros para hacerlo mas ligero y asi poder trabajar con datos
#' mas extensos, en vez de la funcion GAM se utilizo la funcion BAM del paquete mgcv, en vez
#' de la función MAXENT del paquete Dismo se utilizo la función maxnet del paquete maxnet
#' creado por los mismos desarrolladores del MAXENT hecho en JAVA, en vez de la función svm
#' del paquete e1071  se utilizo la función ksvm del paquete kernlab.
#' The LightSSDM package provides five categories of functions (same as SSDM package)(that you can find in details
#' below): Data preparation, Modelling main functions, Model main methods, Model
#' classes, and Miscellaneous.
#'
#' The LightSSDM package provides five categories of functions (same as SSDM package)
#' (that you can find in details below): Data preparation, Modelling main functions,
#' Model main methods, Model classes, and Miscellaneous.
#' \strong{(Esp)} El paquete LightSSDM provee cinco catgorias de funciones (Igual que el paquete SSDM) (pueden encontrarlos en la seccion
#' de abajo): Preparación de datos, funciones principales de modelamiento, métodos principales de modelo,
#' clases de modelo, y varios.
#'
#' @section Data preparation: \describe{\item{\code{\link{load_occ}}}{Load
#'   occurrence data} \item{\code{\link{load_var}}}{Load environmental variables}}
#' @section \strong{(Esp)} Preparación de datos: \describe{\item{\code{\link{load_occ}}}{AbrirLoad
#'   datos de ocurrencia} \item{\code{\link{load_var}}}{Abrir variables ambientales}}
#'
#' @section Modelling main  functions: \describe{
#'   \item{\code{\link{modelling}}}{Build an SDM using a single algorithm}
#'   \item{\code{\link{ensemble_modelling}}}{Build an SDM that assembles
#'   multiple algorithms} \item{\code{\link{stack_modelling}}}{Build an SSDMs
#'   that assembles multiple algorithms and species}}
#' @section \strong{(Esp)} Funciones principales de modelamiento: \describe{
#'   \item{\code{\link{modelling}}}{Construyr un SDM de un solo algoritmo}
#'   \item{\code{\link{ensemble_modelling}}}{Construye un SDM que ensambla
#'   multiples algoritmos} \item{\code{\link{stack_modelling}}}{Construye un SSDMs
#'   que ensambla multiples algoritmos y especies}}
#'
#' @section Model main methods: \describe{
#'   \item{\code{\link{ensemble,Algorithm.SDM-method}}}{Build an ensemble SDM}
#'   \item{\code{\link{stacking,Ensemble.SDM-method}}}{Build an SSDM}
#'   \item{\code{\link{update,Stacked.SDM-method}}}{Update a previous SSDM with
#'   new occurrence data}}
#' @section Métodos principales de modelo: \describe{
#'   \item{\code{\link{ensemble,Algorithm.SDM-method}}}{Construye un SDM ensamblado}
#'   \item{\code{\link{stacking,Ensemble.SDM-method}}}{Construye un SSDM}
#'   \item{\code{\link{update,Stacked.SDM-method}}}{Actualiza un SSDM previó con
#'   nuevos datos de ocurrencias}}
#'
#' @section Model classes: \describe{
#'   \item{\code{\linkS4class{Algorithm.SDM}}}{S4 class to represent SDMs}
#'   \item{\code{\linkS4class{Ensemble.SDM}}}{S4 class to represent ensemble
#'   SDMs} \item{\code{\linkS4class{Stacked.SDM}}}{S4 class to represent SSDMs}}
#' @section  Clases de modelo: \describe{
#'   \item{\code{\linkS4class{Algorithm.SDM}}}{Clase S4 que representa SDMs}
#'   \item{\code{\linkS4class{Ensemble.SDM}}}{Clase S4 que representa un SDM ensamblado}
#'   \item{\code{\linkS4class{Stacked.SDM}}}{Clase S4 que representa SSDMs}}
#'
#' @section Miscellaneous: \describe{ \item{\code{\link{gui}}}{User-friendly
#'   interface for SSDM package} \item{\code{\link{plot.model}}}{Plot SDMs}
#'   \item{\code{\link{save.model}}}{Save SDMs}
#'   \item{\code{\link{load.model}}}{Load SDMs}}
#' @section Varios: \describe{ \item{\code{\link{gui}}}{Interfaz gráfica de usuario
#'   para el paquete LightSSDM} \item{\code{\link{plot.model}}}{Grafica SDMs}
#'   \item{\code{\link{save.model}}}{Guarda SDMs}
#'   \item{\code{\link{load.model}}}{Abre SDMs}}
#'
#' @docType package
#' @name SSDM
NULL
#> NULL
