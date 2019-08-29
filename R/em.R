
#' Expectation maximization algorithm to impute missing GxE cells
#'
#' This function will impute missing cells in the GxE table using an expectation maximization algorithm. 
#'
#' @param Y matrix containing numeric values of cell means with genotypes on rows and environments on columns
#' @param model character vector of length 1. bilinear model to be fit. Arguments can be "AMMI", "GGE", "SREG", "EGE", "GREG". "GGE" and "SREG" are equivalent, as are "EGE" and "GREG".
#' @param tol scalar convergence tolerance threshold, defined as the sum of the absolute value of cell mean differences from iteration i and i-1 scaled by the standard deviation of the values in Y.
#' @param maxiter  integer. Maximum number of iterations.
#' @param k  number of PC to use for imputation. Default is NULL, k will be determined from the imputed data using the parametric bootstrap test. 
#' @param fast logical or integer. If false or 0, k will be deterined at each iteration (slow). If fast is non-zero, k will be estimated each iteration <= max(2, fast), and the last value of k will be used for remaining iterations. . 
#' @param Ytrue  Same as Y but with known, non-mising values. This allows the user to evaluate the accuracy of imputation.
#' @param plotMSE logical. Should the mean square error (MSE) be plotted?.
#' @param verbose logical. Should details e printed?
#' @return Matrix with missing cells replaced by imputed values. 
#' @details
#' 
#' Missing values in the table of genotypes and environments are imputed using an expectation maximization algorithm. The algorithm exits and returns the imputed matrix once a tolerance threshold or maximum number of iterations is reached. This function is generally meant to be used by bilinear when missing cells are found, but the user can also use it to determine the imputation accuracy by providing the true values to 'Ytrue'. 
#' 
#' If 'k' is set to an integer, then this number of PCs will be used for imputation. Otherwise, 'k' will be determined from the model fit using the 'test' argument provided to bilinear. 
#' 
#' If 'fast' is set to TRUE, then the test will only be done for the first 2 iterations. If an integer is provided to 'fast', 'k' will be determined for the first 'fast' iterations.
#' 
#' If a complete matrix of true values is provided, the algorithm will calculate the mean square error. Additionally, if plot MSE is set to true, the MSE of each iteration will be plotted as the algorithm proceeds
#' 
#' If 'verbose' is true, details will be printed to stdout.
#' 
#' @examples
#' data(soy)
#' nMiss <- 10 
#' Ytrue <- soyMeanMat
#' Y <- soyMeanMat
#' Y[sample(1:prod(dim(Y)), nMiss)] <- NA
#'
#' em(Y, model = "AMMI", tol = 1e-5, k = 1, maxiter = 20, Ytrue = Ytrue, plotMSE = TRUE)
#' em(Y, model = "AMMI", tol = 1e-5, k = 2, maxiter = 20, Ytrue = Ytrue, plotMSE = TRUE)
#' em(Y, model = "AMMI", tol = 1e-5, fast = FALSE, maxiter = 20, Ytrue = Ytrue, plotMSE = TRUE)
#' em(Y, model = "AMMI", tol = 1e-5, fast = 2, maxiter = 20, Ytrue = Ytrue, plotMSE = TRUE)
#' @export

em <- function(Y, model, tol = 1e-4, maxiter = 100, k = NULL, fast = TRUE, Ytrue = NULL, plotMSE = FALSE, verbose = FALSE, ...){
	if(!is.null(Ytrue)){
		if (!identical(dim(Y), dim(Ytrue))) stop("To evaluate imputation accuracy using mean square error, please provide a Ytrue of same dimensions as Y.")
		if (!{all(rownames(Ytrue) == rownames(Y)) & all(colnames(Ytrue) == colnames(Y))}) stop("row and col names of Ytrue must match those of Y!")
	}

	mu <- c("(Intercept)" = mean(Y, na.rm = TRUE))
	Eeffect <- colMeans(Y, na.rm = TRUE) - mu
	Geffect <- rowMeans(Y, na.rm = TRUE) - mu

	whichMiss <- which(is.na(Y), arr.ind = TRUE)
	Yimp <- Y
	Yimp[whichMiss] <- mu + Eeffect[whichMiss[, "col"]] + Geffect[whichMiss[, "row"]]

	calcMSE <- !is.null(Ytrue)
	if (calcMSE) {
		Ytrue <- Ytrue[whichMiss]
		mse <- NULL
	}

	sdY <- sd(c(Y), na.rm = TRUE)
	sumdiffscaled <- 1e5
	if (verbose) { if (!is.null(k)) cat("Using", k, "PCs for imputation...\n") else if (fast) cat("Determining number of PCs for imputation...\n") else cat("Determining number of PCs for each iteration...\n") }

	i <- 0
	if (is.null(k) & fast) {
		maxi <- max(2, fast)
		while(i < maxi){
			i <- i + 1
			if (verbose) cat("iteration", i, "\n")
			fitimp <- bilinear(Yimp, verbose = FALSE, returnDataFrame = FALSE, Bonferroni = FALSE, ...)
			k <- fitimp$sigPC
			est <- fitimp$A + fitimp$Theta
			sumdiffscaled <- sum(abs(Yimp[whichMiss] - est[whichMiss])) / sdY
			Yimp[whichMiss] <- est[whichMiss]
			if (calcMSE){
				mse <- c(mse, sum((Yimp[whichMiss] - Ytrue)^2))
			}
		}
		if (k < 1) {
			cat("No significant PCs for imputation! Continuing with 1 PC for imputation...\n")
			k <- 1
		}
		if (verbose) cat("Using", k, "PCs for imputation...\n")
	}

	while(sumdiffscaled > tol & i < maxiter){
		i <- i + 1
		if (fast) {
			fitimp <- bdecomp(Yimp, model = model)
			fitimp[["Theta"]] <- getTheta(k, fitimp$Edecomp)
		} else {
			if (verbose) cat("iteration:", i, "\n")
			fitimp <- bilinear(Yimp, verbose = FALSE, returnDataFrame = FALSE, Bonferroni = FALSE, ...)	
		}
		est <- fitimp$A + fitimp$Theta
		sumdiffscaled <- sum(abs(Yimp[whichMiss] - est[whichMiss])) / sdY
		Yimp[whichMiss] <- est[whichMiss]
		if (calcMSE){
			mse <- c(mse, sum((Yimp[whichMiss] - Ytrue)^2))
		}
		if (verbose & sumdiffscaled < tol) {
			cat("Imputation algorithm converged! returning imputed matrix.\n")
		} else if (i == maxiter) {
			# cat("Imputation algorithm did not converge, maximum iterations reached.\n")
			warning("Imputation algorithm did not converge, maximum iterations reached!")
		}
	}
	if (verbose) cat("total number of iterations:", i, "\n")
	if (plotMSE) plot(1:length(mse), mse, xlab = "iteration", ylab = "Mean Square Error")
	if (calcMSE) return(list(Yimp = Yimp, MSE = mse)) else return(Yimp)
}


