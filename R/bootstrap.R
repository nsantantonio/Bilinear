#' bootstrap function
#'
#' This function perfoms a bootstrap test on one PC at a time inside bilinear(), for ease of parallelization. Not meant to be used by user.
#'
#' @param K integer, kth eigenvalue
#' @param D integer, Maximum right dimensions 
#' @param Dtilde integer. Maximum right dimensions 
#' @param Edecomp svd() object from singular value decomposition of the table of residuals (thetaPlusR)
#' @param M integer. maximum number of PCs.
#' @param B integer. Number of bootstraps to be performed.
#' @param model bilinear model to be fit. Arguments can be "AMMI", "GGE", "SREG", "EGE", "GREG". "GGE" and "SREG" are equivalent, as are "EGE" and "GREG".
#' @param bootMethod character. Method for bootstrap sampling. Can be "full" or "simple", default is "full".
#' @param Theta_k Theta the reduced dimension table can be provided directly. useful for ??
#' @param verbose logical. Should details be printed?
#' @param ... Additional arguments.
#' @return p-values for 1 to K principal components
#' @details
#' 
#' Performs bootstrap test of significance by returning the proportion of time the bootstrap test statistic exceeds the 
#' statistic calculated from the Kth principal component. See Forkman and Piepho (2014) Biometrics, 70(3) for specifics.
#'
#' @keywords parametric bootstrap
#' @importFrom stats rnorm
#' @export
bootstrap <- function(K, D, Dtilde, Edecomp, M, B, model, bootMethod = "full", Theta_k = NULL, verbose = FALSE, ...){
	if (!bootMethod %in% c("full", "simple")) {stop("Please specify the parametric bootstrap method as 'full' or 'simple'.")}
	
	if (verbose & K == 0) cat("Using", bootMethod, "parametric bootstrap method\n")

	Lambda <- Edecomp$d
	nu <- D * Dtilde
	
	I <- nrow(Edecomp$u)
	J <- nrow(Edecomp$v)

	LLt_I <- diag(I) - (1 / I) * matrix(1, I, I) 
	LLt_J <- diag(J) - (1 / J) * matrix(1, J, J) 

	Kplus1 <- K + 1
	Lambda_k <- Lambda[Kplus1:(M)]	

	S <- crossprod(Lambda_k)
	T <- c(Lambda[Kplus1]^2 / S)

	if (bootMethod == "full"){
		sigmasq_k <- 1/nu * S
		if (is.null(Theta_k)){
			if (K == 0){
				Theta_K <- matrix(0, I, J)
			} else if (K == 1){
				Theta_K <- Lambda[1] * tcrossprod(Edecomp$u[,1], Edecomp$v[,1])
			} else {
				Theta_K <- Edecomp$u[, 1:K] %*% tcrossprod(diag(Lambda[1:K]), Edecomp$v[, 1:K])	
			}
		} else {
			Theta_K <- Theta_k[[paste0("PC", K)]]
		}
	}

	T_b <- matrix(NA,nrow = B, ncol = 1)
	b <- 1

	while(b <= B){
		if (bootMethod == "full"){
			Eb <- Theta_K + matrix(rnorm(I * J, sd = sqrt(sigmasq_k)), I, J)
			if (model == "AMMI"){
				Ehatb <- LLt_I %*% Eb %*% LLt_J 
			} else if (model %in% c("GGE","SREG")){
				Ehatb <- LLt_I %*% Eb 
			} else if (model %in% c("EGE","GREG")){
				Ehatb <- Eb %*% LLt_J 
			} else {stop("unspecified model")}
		} else {
			Ehatb<-matrix(rnorm((D - K) * (Dtilde - K)),  nrow = (D - K), ncol = (Dtilde - K))
		}
		
		Lambda_b <- svd(Ehatb)$d
		
		if (bootMethod == "full"){
			T_b[b, 1] <- Lambda_b[Kplus1]^2 / crossprod(Lambda_b[Kplus1:(M)]) 
		} else {
			T_b[b, 1] <- Lambda_b[1]^2 / crossprod(Lambda_b) 
		}
		b <- b + 1
	}
	
	pvalue<-sum(T_b > T)/B
	return(pvalue)
}