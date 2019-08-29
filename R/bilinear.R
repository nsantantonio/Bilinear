#' bilinear function
#'
#' This function will fit bilinear models, such as AMMI, GGE / SREG, EGE / GREG, designed for multienvironment trials from a plant breeding/agronomy perspective, but can be used for fitting any two factor interaction models for dimension reduction of the interaction table. A bootstrap test and an F_R test are implemented for testing the number of significant dimensions.
#'
#' @param x matrix or data.frame containing numeric values of cell means with genotypes on rows and environments on columns (e.g. ith genotype within jth environment). 
#'			Alternatively a data.frame can be provided in long format, with factor and value variable names passed to 'G', 'E' and 'y' respectively
#' @param G character. Name of genotype variable in data.frame 'x' if 'x' is a data.frame in long format. Optionally, character or factor vector of genotype names (optional). If NULL (default) a matrix of cell means must be provided to 'x' argument
#' @param E character. Name of environment variable in data.frame 'x' if 'x' is a data.frame in long format. Optionally, character or factor vector of environment names (optional). If NULL (default) a matrix of cell means must be provided to 'x' argument
#' @param y character. Name of phenotype variable in data.frame 'x' if 'x' is a data.frame in long format. Optionally, numeric vector of phenotype values (optional). If NULL (default) a matrix of cell means must be provided to 'x' argument
#' @param block character. Optional, for RCBD only. Name of block variable in data.frame 'x' if 'x' is a data.frame in long format and study is an RCBD. Optionally, character or factor vector of block names. 
#' @param model bilinear model to be fit. Arguments can be "AMMI", "GGE", "SREG", "EGE", "GREG". "GGE" and "SREG" are equivalent, as are "EGE" and "GREG".
#' @param errorMeanSqDfReps numeric vector containing the error mean square, degrees of freedom, and the number of replications, in that order.
#' @param f numeric value in (0, 1) for scaling exponent on eigenvalues for weighting the genotype scores. Environment scores are weighted 1 - f. Default is 0.5 
#' @param test character. Test for significant PCs. Argument can be either "Ftest" or "bootstrap". Default is "bootstrap".
#' @param imputePC integer. number of PCs for imuptation of missing cells. must be an integer between 1 and min(genotypes, environments) - 2. Default is 1. Alternatively, the argument "sig" can be provided to use the number of significant dimensions from the initial fit (assuming GE_{ij} = 0 for missing cells). THIS HAS NOT YET BEEN IMPLEMENTED!
#' @param alpha pvalue cutoff threshold for significance. Default is 0.05.
#' @param B integer. Number of bootstrap samples for bootstrap test. 1e+06 or greater is recommended, but default is 1e+05 to limit first use computation time.
#' @param nCore integer. Number of cpu cores to be used for parallelization of computation using foreach() and the doMC package. Not tested outside a unix environment. Default is 1.
#' @param Bonferroni logical. Use a bonferroni correction for multiple testing. Default is TRUE
#' @param returnDataFrame logical. Should a data frame of the resulting model fit be returned? Default is TRUE
#' @param override3col logical. Overrides the 3 column numeric warning when a matrix of 3 coulmns is provided. Most users will not need to turn this off, unless the number of environments is 3 and a table of cell means are provided. Default is FALSE.
#' @param verbose logical. When TRUE, details from nalysis are printed.
#' @return The bilinear fit object and tests for significant dimensions.
#' @details
#' 
#' Bilinear models attempt to separate interaction signal from noise in a two factor interaction model by dimension reduction of the table of residuals from the additive model. Dimension reduction is accomplished by singular value decomposition of the residual table, keeping the first k significant dimensions and partitioning the remaining dimensions to the error. 
#' 
#' The interaction model can be formulated for the genotype by environmental interaction as the following: 
#' 
#' y = mu + G_i + E_j + G*E_ij + e_ijk where mu is the grand mean, G_i is the effect of the ith genotype, E_j is the effect of the ith genotype, G*E_ij is the genotype by environment interaction term and e is the experimental error. The G*E_ij term can be thought of as a matrix of cell means (minus the main effects), and can be decomposed as sum(u_ik * d_k * v_jk) + r_ij for $k = 1,...,K$, using singular value decomposition where u_ik is the kth left singular value of the ith genotype, v_jk is the kth right singular value of the jth genotype and d_k is the kth singular value. r is sum(u_ik * d_k * v_jk) for k = K+1, ..., M where M is min(levels(G), levels(E)) - c for c equals 1 or 2 for GGE/SREG/EGE/GREG and AMMI respectively. The genotype and environmental scores are calculated as u_ik * (d_k)^f and v_jk * (d_k)^(1-f) respectively. 
#' 
#' "AMMI" specifies Additive Main Effects Multiplicative Interaction, where $y = mu + G_i + E_j + sum_{k = 1}^K u_{ik} * d_k * v_{jk}) + e$ for $k = 1,...,K$ significant multiplicative terms. 
#' "GGE" or "SREG" specifies a Sites REGression or Genotype + Genotype by Environment interaction, where $y = mu + E_j + sum_{k = 1}^K u_{ik} * d_k * v_{jk}) + e$ for $k = 1,...,K$ significant multiplicative terms. 
#' "EGE" or "GREG" specifies Genotype REGression or Environment + Genotype by Environment interaction,,$y = mu + G_i + sum_{k = 1}^K u_{ik} * d_k * v_{jk}) + e$ for $k = 1,...,K$ significant multiplicative terms. 
#'
#' For specifics see Zobel, R. W., Wright, M. J., & Gauch, H. G. (1988). Statistical analysis of a yield trial. Agronomy Journal, 80(3), 388-393.
#' 
#' Note: either the row and column variable names can specified to 'G' and 'E' respectively and the value variable names to 'y', or by column variable names provided as arguments to 'G', 'E', and 'y'
#' If data.frame already contains variables 'G', 'E' and 'y' ('block' optional) model will be fit accordingly. Any character class in the data.frame is converted to factor.
#' If data is provided in long format and study is replicated within environment, cell means will be used, and a warning will diplay if the data is unbalanced.
#'
#' @examples
#' data(soy)
#' print(soyMeansMat)
#' 
#' AMMIfit <- bilinear(x = soyMeanMat)
#' AMMIfit[c("mu", "Eeffect", "Geffect")] 
#' AMMIfit["DF"] 
#' AMMIfit[c("pvalue", "sigPC")] # Note that the first two terms are considered significant
#' AMMIfit["varcomp"]
#' AMMIfit["svdE"]
#' 
#' # Same soy data in long format
#' # as long as dataframe or matrix column names match c("G", "E", "y") 
#' # exactly (order doesnt matter), the function will assume that the data is in long format
#' print(head(soyMeanDf))
#' AMMIfit <- bilinear(x = soyMeanDf)
#' 
#' # if column names DO NOT match c("G", "E", "y"), they must be supplied to the 'G', 'E', and 'y'
#' # exactly (order doesnt matter), the function will assume that the data is in long format
#' print(head(soyMeanDf))
#' AMMIfit <- bilinear(x = soyMeanDf)
#' # multiple CPUs for bootstrap test, not run
#' # AMMIfit <- bilinear(x = soyMeanMat, model = "AMMI", nCore = 2)
#' 
#' # Sequential F test for significance using cell means. This test is known to be too liberal and will likely result in type I errors
#' AMMIfit <- bilinear(soyMeanMat, test = "Ftest")
#'
#' # F_R test for significance using plot level data in long format to estimate error variance for F test
#' print(head(soy))
#' AMMIfit <- bilinear(soy, test = "Ftest")
#' 
#' # F_R test with estimated error variance. Useful for Augmented designs where the error estimate is estimated from replicated checks
#' # Estimate error variance first 
#' soyfit <- lm(as.formula(y ~ E + G + G:E + block:E), data = soy)
#' errorMeanSq <- summary(soyfit)$sigma^2
#' errorDf <- anova(soyfit)["Residuals", "Df"]
#' reps <- unique(table(soy$G, soy$E))
#' 
#' print(c(errorMeanSq, errorDf, reps))
#' # F_R test with externally estimated error variance
#' soyF2 <- bilinear(soyMeanMat, errorMeanSqDfReps = c(errorMeanSq, errorDf, reps), test = "Ftest")
#' 
#' #reformat data
#' 
#' # as long as dataframe or matrix column names match c("G", "E", "y") 
#' # exactly (order doesnt matter), the function will assume that the data is in long format
#' dfY <- soyMeanDf
#' names(dfY) <- c("G","E","y") # name dfY c("G", "E", "y") 
#' print(head(dfY, 20))
#' bilinear(x = dfY, model = "AMMI", B = 10000, nCore = 2)
#' bilinear(x = soy[c("G","E","y")], model = "AMMI", B = 10000, nCore = 2)
#' 
#' 
#' matY <- as.matrix(dfY)
#' print(matY[1:20,])
#' bilinear(x = matY, model = "AMMI", B = 10000, nCore = 2) 
#' 
#' # if the columns do not match as above, the names of the variables must be
#' # passed to the 'G', 'E', and 'y' arguments as character strings
#' names(dfY) <- c("geno", "env", "trait") #change names dfY
#' print(head(dfY, 20))
#' bilinear(x = dfY, G = "geno", E = "env", y = "trait", model = "AMMI", B = 10000, nCore = 2)
#' data can be also supplied directly as vectors to the 'G', 'E' and 'y' arguments
#' Gvec <- dfY$geno
#' Evec <- dfY$env
#' yvec <- dfY$trait
#' bilinear(G = Gvec, E = Evec, y = yvec, model = "AMMI", alpha = 0.05, B = 10000, nCore = 2)
#'
#' # Impute missing data!
#' missMat <- soyMeanMat 
#' missMat[sample(1:prod(dim(missMat)), 10)] <- NA 
#' print(missMat)
#' bilinear(x = matY, model = "AMMI", B = 10000, nCore = 2) 
#'
#' missMat <- soyMeanMat 
#' missMat[sample(1:prod(dim(missMat)), 10)] <- NA 
#' print(missMat) # 10 missing records
#' missFit <- bilinear(x = missMat, model = "AMMI", B = 10000) 
#' 
#' library(Bilinear)
#' missMeanDf <- soyMeanDf[!{soyMeanDf$E %in% unique(soyMeanDf$E)[4:6] & soyMeanDf$G %in% unique(soyMeanDf$G)[3:4]}, ]
#' nrow(missMeanDf) # 10 missing records
#' missFit2 <- bilinear(x = missMeanDf, model = "AMMI", B = 10000, maxiter = 200, Ytrue = soyMeanMat, plotMSE = TRUE, tol = 1e-4) 
#' @keywords AMMI, GGE, parametric bootstrap
#' @export
bilinear <- function(x = NULL, G = NULL, E = NULL, y = NULL, block = NULL, model = "AMMI", errorMeanSqDfReps = NULL, f=0.5, test = "bootstrap", imputePC = "sig", alpha = 0.05, B = 1e+04, nCore = 1, Bonferroni = FALSE, returnDataFrame = TRUE, override3col = FALSE, verbose = TRUE, ...){

	if (!.Platform$OS.type %in% "unix" & nCore > 1){
		nCore <- 1
		warning("Multicore function only implemented for unix based systems (including macOS and Linux). Continuing with 1 core...\n")
	}

	listOfPackages <- NULL
	if (nCore > 1) listOfPackages <- c(listOfPackages, "foreach")
	
	new_packages <- listOfPackages[!(listOfPackages %in% installed.packages()[, "Package"])]

	if (length(new_packages)){
		cat("To use 'nCore' > 1, please install the 'foreach' package, along your favorite parallel backend that uses 'foreach' (e.g. doMC), load the library and set the number of desired cores (e.g. registerDoMC(nCores)) \n Continuing with nCore = 1\n")	
		nCore <- 1
	}
	
	if (length(listOfPackages)) lapply(listOfPackages, require, character.only = TRUE)

	usrContr <- options("contrasts")[[1]]
	options(contrasts = c("contr.sum", "contr.sum"))
	
	meltName <- function(x, G, E, vName) {
		melted <- data.frame(rep(rownames(x), ncol(x)), rep(colnames(x), each = nrow(x)), c(x))
		names(melted) = c(G, E, vName)
		melted
	}

	renameRank <- function(x){
		x <- t(x)
		nameMat <- matrix(NA, nrow(x), ncol(x), dimnames = list(rownames(x), paste0("rank", 1:ncol(x))))
		for(i in 1:nrow(x)){
			nameMat[i,] <- colnames(x)[order(x[i,])]
		}
		return(nameMat)
	}

	extractPC<-function(k,UDV){
		PC <- c(tcrossprod(UDV$u[,k], UDV$v[,k])) * UDV$d[k]
		return(PC)
	}

	stars <- function(pval){
		ifelse(is.na(pval), "   ", ifelse(pval < 0.001, "***", ifelse(pval < 0.01, " **", ifelse(pval < 0.05, "  *", "   "))))
	}

	anyNullGEy <- is.null(G) | is.null(E) | is.null(y)
	isGEyvec <- any(c(length(G), length(E), length(y)) > 1)
	dataReformatted <- dataFormat(x = x, G = G, E = E, y = y, block = block, alpha = alpha, anyNullGEy = anyNullGEy, override3col = override3col)
	
	Y <- dataReformatted[["Y"]]
	isUnRep <- dataReformatted[["isUnRep"]]
	if (anyNullGEy | isGEyvec) {G <- "G"; E <- "E"; y <- "y"}

	blockSig <- if (!is.null(dataReformatted[["blockSig"]])){if (dataReformatted[["blockSig"]]) TRUE else FALSE} else FALSE

	if (is.null(block) & "block" %in% names(dataReformatted[["DF"]])) block <- "block"
	
	if (!is.null(errorMeanSqDfReps)){
		if (length(errorMeanSqDfReps) != 3) stop("Please provide a vector or list with error variance, error degrees of freedom and number of replications per genotype to the 'errorMeanSqDfReps' argument\n")
		sigmasq <- errorMeanSqDfReps
	} else {
		sigmasq <- dataReformatted[["sigmasq"]]
	}

	if (verbose) {
		cat("Displaying first (up to) 10 rows and columns of two way table of genotype within environment means:\n")
		print(Y[1:min(nrow(Y), 10), 1:min(ncol(Y), 10)])
		cat("\n")
	}

	mu <- c("(Intercept)" = mean(Y, na.rm = TRUE))
	Eeffect <- colMeans(Y, na.rm = TRUE) - mu
	Geffect <- rowMeans(Y, na.rm = TRUE) - mu

	if (verbose) cat("Data consists of", nrow(Y), "Genotypes evaluated in", ncol(Y), "Environments\n")
	I <- nrow(Y) 
	J <- ncol(Y) 

	if (any(is.na(Y))){
		if (verbose)	cat("Missing cells found in the genotype environment table! Imputing missing cells with EM algorithm...\n")
		imp <- em(Y, model = model, verbose = verbose, B = B, ...)
		Y <- if (is.list(imp)) imp$Yimp else imp
	}
	
	decomp <- bdecomp(Y, model = model)
	for (i in names(decomp)) assign(i, decomp[[i]])

	Theta_k <- lapply(0:KmaxPlusOne, getTheta, Edecomp = Edecomp)
	names(Theta_k) <- paste0("PC", 0:KmaxPlusOne)

	sigmasqGxE_k <- c(PC0 = 0)
	sigmasqR_k <- c(PC0 = 1/nu * crossprod(Lambda[1:M]))
	for (k in 1:KmaxPlusOne){
		sigmasqGxE_k[paste0("PC", k)] <- 1/nu * crossprod(Lambda[1:k])
		sigmasqR_k[paste0("PC", k)] <- 1/nu * crossprod(Lambda[(k+1):M])
	}	

	Theta_k <- lapply(Theta_k, function(x){rownames(x) <- rownames(Y); colnames(x) <- colnames(Y); return(x)})
	nominal_k <- lapply(Theta_k, function(x) x + mu + Geffect )

	r <- mean(dataReformatted[["repPerG"]][dataReformatted[["repPerG"]] != 0])
	# r <- mean(c(dataReformatted[["repPerG"]]))
	degfGxE <- c(I + J - (2 * 1:KmaxPlusOne))
	if (model == "AMMI") degfGxE <- degfGxE - 1
	degfR <- nu - cumsum(degfGxE)
 	
 	if (test == "Ftest"){
	 	if (r == 1 & is.null(errorMeanSqDfReps)){
	 		warning("The trial is unreplicated and no error variance was provided, so a sequential F test will be used to determine the number of significant terms. This method is known to be too liberal, and it is suggested that you use the bootstrap method to test for significant terms with the argument test = 'bootstrap'.\n")
	 		Ftest <- (Lambda[1:KmaxPlusOne]^2 / degfGxE)  / ((sigmasqR_k[2:M]) * nu / degfR) 
	 		pvalue <- pf(Ftest, degfGxE, degfR, lower.tail = FALSE)
	 	} else if (r > 1 | !is.null(errorMeanSqDfReps)){
			if (verbose) cat("Using F_R test method\n")
			df1 <- (D - 0:KmaxPlusOne) * (Dtilde - 0:KmaxPlusOne)
			df2 <- sigmasq[[2]]
			reps <- if (r > 1) r else errorMeanSqDfReps[[3]]
			
			Ftest <- c(sigmasqR_k[1:M] * nu) / (df1 * (sigmasq[1] / reps))
			pvalue <- pf(Ftest, df1, df2, lower.tail = FALSE)
		} else {
			stop("oops! something wrong happened >..< please post an issue at https://github.com/nsantantonio/Bilinear/issues ")
		}
	} else if (test == "bootstrap" ){
		if (nCore > 1 & Kmax > 0){
			pvalue <- foreach(K = 0:Kmax, .combine = c) %dopar% bootstrap(K, D = D, Dtilde = Dtilde, M = M, Edecomp = Edecomp, B = B, model = model, verbose = verbose, ...)
		} else {
			pvalue <- NULL
			for (K in 0:Kmax){
				pvalue <- c(pvalue, bootstrap(K, D = D, Dtilde = Dtilde, Edecomp = Edecomp, M = M, B = B, model = model, verbose = verbose, ...))
			}
		}
	} else {
		stop("Please provide either 'bootstrap' or 'Ftest' to test argument.")
	}
	if (Bonferroni){
		if (test == "bootstrap" | (isUnRep & is.null(errorMeanSqDfReps))) alpha <- alpha / 1:KmaxPlusOne else alpha <- alpha / 1:M
	}
	names(pvalue) <- paste0("term", 1:length(pvalue))

	sig <- pvalue < alpha
	Kstar <- if (all(sig)) length(sig) else min(which(!sig)) - 1

	testseq <- if (Bonferroni) " tested sequentially" else ""
	if (verbose) cat(paste0("Number of significant multiplicative terms", testseq, ":"), Kstar, "\n")

	if (nCore > 1){
		PCs <- foreach(k = 1:M, .combine = cbind) %dopar% extractPC(k, UDV = Edecomp)
	} else {
		PCs <- matrix(NA, nrow = I * J, ncol = M)
		for (k in 1:M){
			PCs[,k] <- extractPC(k, UDV = Edecomp)
		}
	}
	colnames(PCs) <- paste0("PC", 1:ncol(PCs))

	Theta <- Theta_k[[paste0("PC", Kstar)]]

	if (Kstar > 0){
		if (Kstar == 1){
			Gscores <- Edecomp$u[, 1, drop = FALSE] * Lambda[1]^(f)
			Escores <- Edecomp$v[, 1, drop = FALSE] * Lambda[1]^(1-f)
		} else {
			Gscores <- Edecomp$u[, 1:Kstar] %*% diag(Lambda[1:Kstar]^(f))
			Escores <- Edecomp$v[, 1:Kstar] %*% diag(Lambda[1:Kstar]^(1-f))
		}
		colnames(Escores) <- colnames(Gscores) <- paste0("PC", 1:Kstar)
	} else {
		Gscores <- matrix(0, I)
		Escores <- matrix(0, J)
		colnames(Escores) <- colnames(Gscores) <- "NoSigGxE"
	}
	
	rownames(Gscores) <- names(Geffect)
	rownames(Escores) <- names(Eeffect)

	scores <- list(Gscores = Gscores, Escores = Escores)		

	rankTables <- lapply(nominal_k[1:(Kstar + 1)], function(x) apply(-x, 2, rank))
	nameRank <- lapply(rankTables, renameRank)	

	Sigmasq <- c(sigmasqGxE = sigmasqGxE_k[[paste0("PC", Kstar)]], sigmasqR = sigmasqR_k[[paste0("PC", Kstar)]])
	if (!is.null(sigmasq)) {Sigmasq["sigmasq"] <- sigmasq[1]}
	degf <- c(degfGxE = sum(degfGxE[1:Kstar]), degfR = degfR[Kstar])
	if (!is.null(sigmasq)) {degf["degfErr"] <- sigmasq[2]}

	rownames(Theta) <- rownames(Y) 
	colnames(Theta) <- colnames(Y) 
	R <- Y - A - Theta

	addEffects <- list(mu = mu, Eeffect = Eeffect, Geffect = Geffect)

	if (returnDataFrame){
		coeff.list <- list(Y = meltName(Y, G = G, E = E, vName = "Y"),
						   A = meltName(A, G = G, E = E, vName = "A"), 
						   ThetaPlusR = meltName(Emat, G = G, E = E, vName = "ThetaPlusR"), 
						   Theta = meltName(Theta, G = G, E = E, vName = "Theta"), 
						   R = meltName(R, G = G, E = E, vName = "R"))
		fitDF <- data.frame(Reduce(function(...) merge(..., by = c(E, G)), coeff.list), PCs)
		output <- list(DF = fitDF)
	} else {
		PClist <- sapply(colnames(PCs), function(x) NULL)
		for(i in 1:ncol(PCs)){
			PClist[[i]] <- matrix(PCs[,i], I, J)
		}
		output <- c(list(Y = Y, A = A, Theta = Theta, R = R), PClist)
	}

	if (isUnRep) nPCANOVA <- KmaxPlusOne else nPCANOVA <- M
	
	blockStatement <- if (blockSig) paste0(" + ", E, ":", block) else ""

	modelStatement <- paste0(y, " ~ ", E, blockStatement, " + ", G)
	if (!isUnRep) modelStatement <- paste0(modelStatement, " + ", G, ":", E)
	fit <- lm(as.formula(modelStatement), data = dataReformatted[["DF"]])

	anovafit <- as.data.frame(anova(fit)) 
	names(anovafit) <- c("Df", "SS", "MS", "testStatistic", "Pvalue")
	anovaDf <- c(degfGxE, tail(degfR, 1))

	SS <- Lambda[1:M]^2
	if (r > 1) SS <- SS * r
	MS <- SS / anovaDf
	PCpvalue <- pvalue
	Tstat <- if (test == "Ftest") Ftest else pvalue * B
	if ({test == "bootstrap" | isUnRep} & is.null(errorMeanSqDfReps)) {
		PCpvalue <- c(PCpvalue, NA)
		Tstat <- c(Tstat, NA)
	}
	names(PCpvalue) <- paste0("PC", 1:length(PCpvalue))
	names(Tstat) <- paste0("PC", 1:length(Tstat))
	anovaPC <- data.frame(Df = anovaDf, SS = SS, MS = MS, testStatistic = Tstat, Pvalue = PCpvalue)
	if (isUnRep){
		ANOVA <- rbind(anovafit[c(E, G),], anovaPC)
	} else if (blockSig) {
		ANOVA <- rbind(anovafit[c(E, G, paste0(E, ":", block)),], anovaPC, anovafit["Residuals", ])
	} else {
		ANOVA <- rbind(anovafit[c(E, G),], anovaPC, anovafit["Residuals", ])
	}

	if (blockSig){
		ANOVA[E, "testStatistic"] <- ANOVA[E, "MS"] / ANOVA[paste0(E, ":", block), "MS"]
		ANOVA[E, "Pvalue"] <- pf(ANOVA[E, "testStatistic"], df1 = ANOVA[E, "Df"], df2 = ANOVA[paste0(E, ":", block), "MS"], lower.tail = FALSE)
	}

	if (!is.null(errorMeanSqDfReps)) {
		rownames(ANOVA)[nrow(ANOVA)] <- paste0("PC", M)
		ANOVA <- rbind(ANOVA, data.frame(Df = errorMeanSqDfReps[[2]], SS = prod(unlist(errorMeanSqDfReps[1:2])), MS = errorMeanSqDfReps[[1]], testStatistic = NA, Pvalue = NA, row.names = "Residuals"))		
		for (i in c(G, E)) {
			ANOVA[i, "testStatistic"] <- ANOVA[i, "MS"] / ANOVA["Residuals", "MS"]
			ANOVA[i, "Pvalue"] <- pf(ANOVA[i, "testStatistic"], df1 = ANOVA[i, "Df"], df2 = ANOVA["Residuals", "MS"], lower.tail = FALSE)
		}
	}

	if (verbose) {
		printANOVA <- ANOVA
		printANOVA[[" "]] <- stars(printANOVA[["Pvalue"]])
		printANOVA[["Pvalue"]] <- sprintf("%e", printANOVA[["Pvalue"]]) 
		printANOVA[printANOVA[["Pvalue"]] %in% "NA", "Pvalue"] <- NA 

		printANOVA[as.numeric(printANOVA[["Pvalue"]]) %in% 0, "Pvalue"] <- paste0("< ", 1 / B)
		printANOVA[is.na(printANOVA[["Pvalue"]]), "Pvalue"] <- ""
		cat("Analysis of Variance Table\n", paste0("Response: ", y), "\n")
		print(printANOVA)
		cat("---\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05\n")
		cat("Number of significant multiplicative terms (tested sequentially): ", Kstar, "\n")
		if (blockSig) cat("NOTE: P-values for additive environment effects are tested with the", paste0(E,":",block), "term.\n")
	}
	options(contrasts = usrContr)
	return(c(model = model, addEffects, list(ANOVA = ANOVA), output, list(pvalue = pvalue, sigPC = Kstar, scores = scores, varcomp = Sigmasq,
	 		 degf = degf, svdE = Edecomp, rankTables = rankTables, nameRank = nameRank)))
}

