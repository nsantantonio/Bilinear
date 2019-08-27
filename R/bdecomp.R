#' Expectation maximization algorithm to impute missing GxE cells
#'
#' This function will impute missing cells in the GxE table using an expectation maximization algorithm. Not intended for direct use by the user.
#'
#' @param Y matrix containing numeric values of cell means with genotypes on rows and environments on columns 
#' @param model bilinear model to be fit. Arguments can be "AMMI", "GGE", "SREG", "EGE", "GREG". "GGE" and "SREG" are equivalent, as are "EGE" and "GREG".
#' @return Returns the bilinear decomposition of the GE matrix 
#' @details
#' 
#' A: Matrix of additive effects (mu + G + E)
#' Emat: Matrix of residuals from the additive mode
#' Edecomp: singlualr value decomposition of Emat
#' Lambda: singular values from 
#' nu: degrees of freedom
#' 
#' @examples
#' data(soy)
#' bdecomp(soyMeansMat, "AMMI")
#' @export

bdecomp <- function(Y, model) {

	I <- nrow(Y) 
	J <- ncol(Y) 
	# if (verbose) cat("Data consists of ", I, " Genotypes evaluated in ", J, " Environments\n")
	
	if (model == "AMMI"){
		LLt_I <- diag(I) - (1 / I) * matrix(1, I, I) 
		LLt_J <- diag(J) - (1 / J) * matrix(1, J, J) 
		M <- min(I - 1, J - 1)
		D <- I - 1
		Dtilde <- J - 1
		Emat <- LLt_I %*% Y %*% LLt_J
	} else if (model %in% c("GGE","SREG")){
		LLt_I <- diag(I) - (1 / I) * matrix(1, I, I) 
		M <- min(I - 1, J)
		D <- I - 1 
		Dtilde <- J
		Emat <- LLt_I %*% Y
	} else if (model %in% c("EGE","GREG")){
		LLt_J <- diag(J) - (1 / J) * matrix(1, J, J) 
		M <- min(I, J - 1)
		D <- I
		Dtilde <- J - 1
		Emat <- Y %*% LLt_J
	} else {
		stop("'model' argument must be either 'AMMI', 'GGE', 'SREG', 'EGE',or 'GREG'!\nNote: 'GGE' and 'SREG' are equivalent, as are 'EGE' and 'GREG'.")
	}

	rownames(Emat) <- rownames(Y) 
	colnames(Emat) <- colnames(Y) 
	nu <- D * Dtilde
	
	Kmax <- M - 2
	KmaxPlusOne <- Kmax + 1
	A <- Y - Emat

	Edecomp <- svd(Emat)
	Lambda <- Edecomp$d

	list(A = A, Emat = Emat, Edecomp = Edecomp, Lambda = Lambda, nu = nu, M = M, D = D, Dtilde = Dtilde, Kmax = Kmax, KmaxPlusOne = KmaxPlusOne)
}
