#' Extract Theta from decomposition of GxE residual matrix from additive, for some k in {0, 1, ..., M}.
#'
#' Calculate dimension reduction of kth degree 
#'
#' @param k number of PCs, k in {0, 1, ..., M}.
#' @param Edecomp decomposition of GxE residual matrix from additive
#' @return Returns the matrix Theta of dimension k
#' @details
#' 
#' Theta: matrix of dimension k, from decomposition of GxE residual 
#' 
#' @examples
#' data(soy)
#' Edecomp <- bdecomp(soyMeansMat, "AMMI")
#' getTheta(Edecomp, 2)
#'
#' @export

getTheta <- function(k, Edecomp){
	I <- nrow(Edecomp$u)
	J <- nrow(Edecomp$v)
	Theta <- if(k == 0){
		matrix(0, I, J)
	} else if(k == 1) {
		Edecomp$d[1] * tcrossprod(Edecomp$u[,1], Edecomp$v[,1])
	} else {
		Edecomp$u[, 1:k] %*% tcrossprod(diag(Edecomp$d[1:k]), Edecomp$v[, 1:k])
	}
	Theta
}