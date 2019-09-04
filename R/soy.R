#' Soybean data from Zobel, Wright and Gauch (1988) and 
#' is subset from the \code{"gauch.soy"} dataset included in the \code{"agridat"} package 
#' containing 7 genotypes evaluated in 10 environments
#' 
#' @name soy
#' @docType data
#'
#' @usage data(soy)
#'
#' @format 
#' soy An object of class \code{"data.frame"} containing plot values. The blocks (or reps) are randomly assigned as explained from the documentation of \code{"gauch.soy"} data from \code{"agridat"}
#' 
#' @keywords soy
#'
#' @references Zobel, R. W., Wright, M. J., & Gauch, H. G. (1988). 
#' Statistical analysis of a yield trial. Agronomy Journal, 80(3), 388-393.
#' (\href{https://dl.sciencesocieties.org/publications/aj/abstracts/80/3/AJ0800030388}{Agronomy Journal})
#'
#' @source \href{https://cran.r-project.org/web/packages/agridat/index.html}{agridat}
#'
#' @examples
#' data(soy)
#' summary(soy)
#' print(soy)
NULL
