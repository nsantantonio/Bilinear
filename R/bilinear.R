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
#' @return The bilinear fit object and tests for significant dimensions.
#' @details
#' 
#' Bilinear models attempt to separate interaction signal from noise in a two factor interaction model by dimension reduction of the table of residuals from the additive model. Dimension reduction is accomplished by singular value decomposition of the residual table, keeping the first k significant dimensions and partitioning the remaining dimensions to the error. 
#' 
#' The interaction model can be formulated for the genotype by environmental interaction as the following: 
#' 
#' y = mu + G_i + E_j + G*E_ij + e_ijk where mu is the grand mean, G_i is the effect of the ith genotype, E_j is the effect of the ith genotype, G*E_ij is the genotype by environment interaction term and e is the experimental error. The G*E_ij term can be thought of as a matrix of cell means (minus the main effects), and can be decomposed as sum(u_ik * d_k * v_jk) + r_ij for $k = 1,...,K$, using singular value decomposition where u_ik is the kth left singular value of the ith genotype, v_jk is the kth right singular value of the jth genotype and d_k is the kth singular value. r is sum(u_ik * d_k * v_jk) for k = K+1, ..., M where M is min(levels(G), levels(E)) - c for c equals 1 or 2 for GGE/SREG/EGE/GREG and AMMI respectively. The genotype and environmental scores are calculated as u_ik * (d_k)^f and v_jk * (d_k)^(1-f) respectively. 
#' 
#' "AMMI" specifies Additive Main Effects Multiplicative Interaction, where $y = \mu + G_i + E_j + \sum_{k = 1}^K u_{ik} * d_k * v_{jk}) + e$ for $k = 1,...,K$ significant multiplicative terms. 
#' "GGE" or "SREG" specifies a Sites REGression or Genotype + Genotype by Environment interaction, where $y = \mu + E_j + \sum_{k = 1}^K u_{ik} * d_k * v_{jk}) + e$ for $k = 1,...,K$ significant multiplicative terms. 
#' "EGE" or "GREG" specifies Genotype REGression or Environment + Genotype by Environment interaction,,$y = \mu + G_i + \sum_{k = 1}^K u_{ik} * d_k * v_{jk}) + e$ for $k = 1,...,K$ significant multiplicative terms. 
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
#' soyF2 <- bilinear(soyShort, errorMeanSqDfReps = c(errorMeanSq, errorDf, reps), test = "Ftest")
#' 
#' #reformat data
#' 
#' # as long as dataframe or matrix column names match c("G", "E", "y") 
#' # exactly (order doesnt matter), the function will assume that the data is in long format
#' names(dfY) <- c("G","E","y") # name dfY c("G", "E", "y") 
#' print(head(dfY, 20))
#' bilinear(x = dfY, model = "AMMI", B = 10000, nCore = 2)
#' bilinear(x = soyLong[c("G","E","y")], model = "AMMI", B = 10000, nCore = 2)
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
# data can be also supplied directly as vectors to the 'G', 'E' and 'y' arguments
#' Gvec <- dfY$geno
#' Evec <- dfY$env
#' yvec <- dfY$trait
#' bilinear(G = Gvec, E = Evec, y = yvec, model = "AMMI", alpha = 0.05, B = 10000, nCore = 2)
#' @keywords AMMI, GGE, parametric bootstrap
#' @export
bilinear <- function(x = NULL, G = NULL, E = NULL, y = NULL, block = NULL, model = "AMMI", errorMeanSqDfReps = NULL, f=0.5, test = "bootstrap", imputePC = "sig", alpha = 0.05, B = 1e+04, nCore = 1, Bonferroni = TRUE, returnDataFrame = TRUE, override3col = FALSE, ...){
# x = Y; f = 0.5; G = NULL; E = NULL; y = NULL; block = NULL; model = "AMMI"; test = "bootstrap"; errorMeanSqDfReps = NULL; alpha = 0.05; B = 10000; nCore = 2; Bonferroni = TRUE; returnDataFrame = TRUE; override3col = TRUE;
	

	if(!.Platform$OS.type %in% "unix" & nCore > 1){
		nCore <- 1
		warning("Multicore function only implemented for unix based systems (including OS X and Linux).
Continuing with 1 core... if analysis takes too long, press Ctrl+c to quit analysis.")
	}

	list_of_packages <- c("reshape2")
	if(nCore > 1) list_of_packages <- c(list_of_packages, "foreach", "doMC")
	
	new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

	if(length(new_packages)){
		cat("Additional Packages required to run:\n",list_of_packages, "\nDo you want to install them?\n")	
		install.packages(new_packages)
	}
	
	lapply(list_of_packages, require, character.only = TRUE)
	if(nCore > 1) registerDoMC(cores = nCore)
	

	# list.of.packages <- c("reshape2","foreach", "doMC")
	# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	# if(length(new.packages)) {install.packages(new.packages)}
	
	# packages.required<-c(TRUE, TRUE, nCore > 1)
	
	# lapply(list.of.packages[packages.required], require, character.only = TRUE)

	# if(nCore > 1) {registerDoMC(cores = nCore)}

	usrContr <- options("contrasts")[[1]]
	options(contrasts = c("contr.sum", "contr.sum"))
	
	meltName <- function(x, vName){
		x <- melt(x)
		names(x) <- c(G, E, vName)
		return(x)
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

	bootstrap <- function(K, D, Dtilde, Edecomp, B, bootMethod = "full", Theta_k = NULL, ...){
		if(!bootMethod %in% c("full", "simple")) {stop("Please specify the parametric bootstrap method as 'full' or 'simple'.")}
		
		if(K == 1) cat("Using", bootMethod, "parametric bootstrap method\n")

		Lambda <- Edecomp$d
		nu <- D * Dtilde
		
		I <- nrow(Edecomp$u)
		J <- nrow(Edecomp$v)

		Kplus1 <- K + 1
		Lambda_k <- Lambda[Kplus1:(M)]	

		S <- crossprod(Lambda_k)
		T <- c(Lambda[Kplus1]^2 / S)

		if(bootMethod == "full"){
			sigmasq_k <- 1/nu * S
			if(is.null(Theta_k)){
				if (K == 0){
					Theta_K <- matrix(0, I, J)
				} else if(K == 1){
					Theta_K <- Lambda[1] * tcrossprod(Edecomp$u[,1], Edecomp$v[,1])
				} else {
					Theta_K <- Edecomp$u[, 1:K] %*% tcrossprod(diag(Lambda[1:K]), Edecomp$v[, 1:K])	
				}
			} else {
				Theta_K <- Theta_k[[paste0("PC", k)]]
			}
		}

		T_b <- matrix(NA,nrow = B, ncol = 1)
		b <- 1

		while(b <= B){
			if(bootMethod == "full"){
				Eb <- Theta_K + matrix(rnorm(I * J, sd = sqrt(sigmasq_k)), I, J)
				if(model == "AMMI"){
					Ehatb <- LLt_I %*% Eb %*% LLt_J 
				} else if(model %in% c("GGE","SREG")){
					Ehatb <- LLt_I %*% Eb 
				} else if(model %in% c("EGE","GREG")){
					Ehatb <- Eb %*% LLt_J 
				} else {stop("unspecified model")}
			} else {
				Ehatb<-matrix(rnorm((D - K) * (Dtilde - K)),  nrow = (D - K), ncol = (Dtilde - K))
			}
			
			Lambda_b <- svd(Ehatb)$d
			
			if(bootMethod == "full"){
				T_b[b, 1] <- Lambda_b[Kplus1]^2 / crossprod(Lambda_b[Kplus1:(M)]) 
			} else {
				T_b[b, 1] <- Lambda_b[1]^2 / crossprod(Lambda_b) 
			}
			
			b <- b + 1
		}
		
		pvalue<-sum(T_b > T)/B
		return(pvalue)
	}

	dataFormat <- function(...){

		warnmessage <- "Please provide a dataframe or matrix for x!
If in wide format, put genotypes as rows and environments
as columns. If in long format, provide the names of the dataframe 
variables that correspond to genotypes, environments and the trait value 
as character strings to the 'G', 'E' and 'y' arguments, respectively.
If the experiment is a randomized complete block design, provide
the name of the blocking variable to 'block'.

Alternatively, you may provide vectors of genotype names, 
environment names and the trait value directly to the 'G', 
'E', and 'y' arguments, respectively. If applicable also 
provide a vector of blocks to the 'block' argument\n"

		warnlackrank <- "At least 3 Genotypes and 3 Environments are necessary for these bilinear models!\n"

		warnunbal <- "WARNING: The experimental design is unbalanced in that not all genotypes were replicated 
the same number of times within each environment.  Means within sites will be
used since replications differ. Care should be taken when drawing inference
from results of unbalanced designs.\n"

		warnMissingCells <- "WARNING: At least one genotype was not observed in at least one environment.
An EM algorithm will be used to estimate the missing cells. 
Care should be taken when interpreting results from designs with 
missing genotype environment combinations.\n"

		blockNotSig <- "Block was not significant at specified alpha level, 
continuing using the pooled error term, or value provided to errorMeanSqDfReps.  
To add a different alpha level for the block effect provide a vector of 
length 2 to the alpha argument with the first element as the significance 
threshold for the number of PCs,and the second as the significance threshold
for the 'block' effect.\n"

		isRCBD <- !is.null(block) | "block" %in% colnames(x)
		anyNullGEy <- is.null(G) | is.null(E) | is.null(y)
		allNullGEy <- is.null(G) & is.null(E) & is.null(y)
		anyMissCells <- NULL

		if (is.null(x)){
			if (anyNullGEy | !is.numeric(y)){
				stop(warnmessage)
			} else {
				DF <- data.frame(G = G, E = E, y = y)
				G <- "G"; E <- "E";	y <- "y"
				if (isRCBD) {
					DF$block <- block
					block <- "block"
				}
			}
		} else {
			isAllNum <- all(sapply(x,class) %in% c("numeric", "integer")) 

			if (ncol(x) == 3 & anyNullGEy & isAllNum){
				is3env <- "yes"
				if(!override3col) is3env <- readline("You have provided a dataframe or matrix with 3 columns, and no variable names to 'G', 'E' and/or 'y'.
							Is the provided dataframe or matrix contain 3 environments with genotypes on the rows?  (y/n) ")
				if(override3col | is3env %in% c("y","yes","YES","Yes","Y")) {
					x <- as.matrix(x)
				} else {
					stop(warnmessage)
				}
			}

			if (all(c("G","E","y") %in% colnames(x))){
				DF <- as.data.frame(x[,c("G", "E", "y")])
				G <- "G"; E <- "E";	y <- "y"
				if(!class(DF[[y]]) %in% "numeric") DF[[y]] <- as.numeric(DF[[y]])
				if(isRCBD) {
					DF$block <- x[,"block"]
					block <- "block"
				}
			} else  if (is.data.frame(x)){
				DF <- x
			} else if (is.matrix(x)){
				if (is.numeric(x) & allNullGEy) {
					if(any(is.na(x))) anyMissCells <- which(is.na(x), arr.ind = TRUE)
					DF <- melt(x)
					names(DF) <- c("G", "E", "y")
					return(list(Y = x, DF = DF, isUnRep = TRUE, sigmasq = NULL, repPerG = 1, anyMissCells = anyMissCells))
				} else if(!allNullGEy){
					DF <- x
				}
			} else {
				stop(warnmessage)
			}			
		}
		

		Elvls <- unique(DF[[E]][!is.na(DF[[E]])])
		Glvls <- unique(DF[[G]][!is.na(DF[[G]])])

		DF[[E]] <- factor(DF[[E]], levels = Elvls)
		DF[[G]] <- factor(DF[[G]], levels = Glvls)

		repGE <- table(DF[[G]][!is.na(DF[[y]])], DF[[E]][!is.na(DF[[y]])])
		nReps <- unique(repGE)
		isUnRep <- length(nReps) == 1 & nReps[1] == 1
		
		I <- nrow(repGE)
		J <- ncol(repGE)

		if(I < 3 | J < 3) {stop(warnlackrank)}
		
		if (length(nReps) > 1){
			cat(warnunbal)
			if (any(repGE == 0)){
				stop("As of now this program cannot handle missing cells in the genotype environment table. 
Please provide a dataset with all genotypes observedin all environments. 
An EM algorithm is planned to be implemented for this purpose and will be released in a later version.")
				 # cat(warnMissingCells)
				# anyMissCells <- which(repGE == 0, arr.ind = TRUE)
				# DFmu <- aggregate(as.formula(paste0(y," ~ ", E, " + ", G)), FUN = mean, data = DF)
				# Y <- acast(DFmu,as.formula(paste0(G, " ~ ", E)), value.var = y)
				# mu <- mean(Y, na.rm = TRUE)
				# Gmu <- rowMeans(Y, na.rm = TRUE) - mu
				# Emu <- colMeans(Y, na.rm = TRUE) - mu
				# GEmu <- Gmu
			}
		} 

		if(isUnRep){
			fit <- lm(as.formula(paste0(y," ~ ", E, " + ", G)), data = DF)
			Y <- acast(DF, as.formula(paste0(G, " ~ ", E)), value.var = y)
		} else {
			if (isRCBD){	
				fit <- lm(as.formula(paste0(y," ~ ", E, " + ", G, " + ", G, ":", E, " + ", block , ":", E)), data = DF)
				pvalBlock <- anova(fit)["E:block","Pr(>F)"]
				blockSig <- TRUE
				if (length(alpha) > 1){
					alphaBlock <- alpha[2] 
					alpha <- alpha[1]
				} else {
					alphaBlock <- alpha
				}
				if (pvalBlock >= alphaBlock){
					cat(blockNotSig)
					blockSig <- FALSE
					fit <- lm(as.formula(paste0(y," ~ ", E, " + ", G, " + ", G, ":", E)), data = DF)
				}	
			} else {
				fit <- lm(as.formula(paste0(y," ~ ", E, " + ", G, " + ", G, ":", E)), data = DF)
			}
			coefs <- dummy.coef(fit)
			intercept <- coefs[[1]]

			Emat <- matrix(rep(coefs[[2]], I), nrow = I, ncol = J, byrow = TRUE)
			Gmat <- matrix(rep(coefs[[3]], J), nrow = I, ncol = J)
			GEmat <- matrix(coefs[[4]], nrow = I, ncol = J, byrow = TRUE)
			
			Y <- intercept + Gmat + Emat + GEmat
			rownames(Y) <- Glvls
			colnames(Y) <- Elvls
		}
		if(isUnRep) sigmasq <- NULL else sigmasq <- summary(fit)$sigma^2
		errorDf <- anova(fit)["Residuals", "Df"]
		returnList <- list(Y = Y, DF = DF, isUnRep = isUnRep)

		if(!isUnRep){
			sigmaList <- c(sigmasq = sigmasq, errorDf = errorDf)
			returnList <- c(returnList, list(sigmasq = sigmaList))
		}
		if(length(nReps) > 1) repPerG <- repGE else repPerG <- nReps
		if(isRCBD) returnList <- c(returnList, blockSig = blockSig)
		returnList <- c(returnList, list(repPerG = repPerG), anyMissCells = anyMissCells)
		return(returnList)
	}

	anyNullGEy <- is.null(G) | is.null(E) | is.null(y)
	isGEyvec <- any(c(length(G), length(E), length(y)) > 1)
	dataReformatted <- dataFormat(x = x, G = G, E = E, y = y, block = block, alpha = alpha, anyNullGEy = anyNullGEy)
	
	Y <- dataReformatted[["Y"]]
	isUnRep <- dataReformatted[["isUnRep"]]
	if(anyNullGEy | isGEyvec) {G <- "G"; E <- "E"; y <- "y"}
	if(is.null(block) & "block" %in% names(dataReformatted[["DF"]])) block <- "block"
	
	if(!is.null(errorMeanSqDfReps)){
		if(length(errorMeanSqDfReps) != 3) stop("Please provide a vector with error variance, 
error degrees of freedom and number of replications per genotype to the 'errorMeanSqDfReps' argument\n")
		sigmasq <- errorMeanSqDfReps
	} else {
		sigmasq <- dataReformatted[["sigmasq"]]
	}


	cat("\nDisplaying first (up to) 10 rows and columns of two way table of genotype within environment means:\n\n")
	print(Y[1:min(nrow(Y), 10), 1:min(ncol(Y), 10)])
	cat("\n")

	mu <- c("(Intercept)" = mean(Y, na.rm = TRUE))
	Eeffect <- colMeans(Y, na.rm = TRUE) - mu
	Geffect <- rowMeans(Y, na.rm = TRUE) - mu

	I <- nrow(Y) 
	J <- ncol(Y) 
	cat("\nData consists of ", I, " Genotypes evaluated in ", J, " Environments\n\n")
	
	# fitBilinear <- (...){

	if(model == "AMMI"){
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
		stop("'model' argument must be either 'AMMI', 'GGE', 'SREG', 'EGE',or 'GREG'!
Note: 'GGE' and 'SREG' are equivalent, as are 'EGE' and 'GREG'.")
	}

	rownames(Emat)<-rownames(Y) 
	colnames(Emat)<-colnames(Y) 
	nu <- D * Dtilde
	
	Kmax <- M - 2
	KmaxPlusOne <- Kmax + 1
	A <- Y - Emat

	Edecomp <- svd(Emat)
	Lambda <- Edecomp$d

	Theta_k <- list(PC0 = matrix(0, I, J))
	sigmasqGxE_k <- c(PC0 = 0)
	sigmasqR_k <- c(PC0 = 1/nu * crossprod(Lambda[1:M]))

	sigmasqGxE_k["PC1"] <- 1/nu * Lambda[1]^2 
	sigmasqR_k["PC1"] <- 1/nu * crossprod(Lambda[2:M]) 
	Theta_k[["PC1"]] <- Lambda[1] * tcrossprod(Edecomp$u[,1], Edecomp$v[,1])

	if(KmaxPlusOne > 1){
		for (k in 2:KmaxPlusOne){
			sigmasqGxE_k[paste0("PC", k)] <- 1/nu * crossprod(Lambda[1:k])
			sigmasqR_k[paste0("PC", k)] <- 1/nu * crossprod(Lambda[(k+1):M])
			Theta_k[[paste0("PC", k)]] <- Edecomp$u[, 1:k] %*% tcrossprod(diag(Lambda[1:k]), Edecomp$v[, 1:k])
		}
	}	
	# sigmasqR_k[paste0("PC", M)] <- 1/nu * Lambda[M]^2
	Theta_k <- lapply(Theta_k, function(x){rownames(x) <- rownames(Y); colnames(x) <- colnames(Y); return(x)})
	nominal_k <- lapply(Theta_k, function(x) x + mu + Geffect )
 
# F tests for significance

	n <- mean(c(dataReformatted[["repPerG"]]))
	degfGxE <- c(I + J - (2 * 1:KmaxPlusOne))
	if(model == "AMMI") degfGxE <- degfGxE - 1
	degfR <- nu - cumsum(degfGxE)
 	
 	if(test == "Ftest"){
	 	if(n == 1 & is.null(errorMeanSqDfReps)){
	 		warning("The trial is unreplicated and no error variance was provided,
so a sequential F test will be used to determine the number of significant terms.
This method is known to be too liberal, and it is suggested that you use the bootstrap
method to test for significant terms with the argument test = 'bootstrap'.\n")
	 		Ftest <- (Lambda[1:KmaxPlusOne]^2 / degfGxE)  / ((sigmasqR_k[2:M]) * nu / degfR) 
	 		pvalue <- pf(Ftest, degfGxE, degfR, lower.tail = FALSE)
	 	} else if(n > 1 | !is.null(errorMeanSqDfReps)){
			cat("Using F_R test method\n")
			df1 <- (D - 0:KmaxPlusOne) * (Dtilde - 0:KmaxPlusOne)
			df2 <- sigmasq[[2]]
			
			if(n > 1) reps <- n else reps <- errorMeanSqDfReps[3]
			
			Ftest <- c(sigmasqR_k[1:M] * nu) / (df1 * (sigmasq[1] / reps))
			pvalue <- pf(Ftest, df1, df2, lower.tail = FALSE)
		} else {stop("oops! something wrong happened >..<")}
	} else if (test == "bootstrap" ){
		if(nCore > 1 & Kmax > 0){
			pvalue <- foreach(K = 0:Kmax, .combine = c) %dopar% bootstrap(K, D = D, Dtilde = Dtilde, Edecomp = Edecomp, B = B, ...)
		} else {
			pvalue <- bootstrap(0, D = D, Dtilde = Dtilde, Edecomp = Edecomp, B = B, ...)
			if(Kmax > 0){
				for (K in 1:Kmax){
					pvalue <- c(pvalue, bootstrap(K, D = D, Dtilde = Dtilde, Edecomp = Edecomp, B = B, ...))
				}
			}
		}
	}
# browser()
	if(Bonferroni){
		if(test == "bootstrap" | (isUnRep & is.null(errorMeanSqDfReps))) alpha <- alpha / 1:KmaxPlusOne else alpha <- alpha / 1:M	
	}
	names(pvalue) <- paste0("term", 1:length(pvalue))

	sig <- pvalue < alpha
	if(all(sig)) Kstar <- length(sig) else Kstar <- min(which(!sig)) - 1

	# cat("P-values of multiplicative terms: \n"); print(pvalue); cat("\n")
	cat("Number of significant multiplicative terms (tested sequentially): ", Kstar, "\n")

	if(nCore > 1){
		PCs <- foreach(k = 1:M, .combine = cbind) %dopar% extractPC(k, UDV = Edecomp)
	} else {
		PCs <- matrix(NA, nrow = I * J, ncol = M)
		for (k in 1:M){
			PCs[,k] <- extractPC(k, UDV = Edecomp)
		}
	}
	colnames(PCs) <- paste0("PC", 1:ncol(PCs))

	Theta <- Theta_k[[paste0("PC", Kstar)]]

	if(Kstar > 0){
		if(Kstar == 1){
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
	if(!is.null(sigmasq)) {Sigmasq["sigmasq"] <- sigmasq[1]}
	degf <- c(degfGxE = sum(degfGxE[1:Kstar]), degfR = degfR[Kstar])
	if(!is.null(sigmasq)) {degf["degfErr"] <- sigmasq[2]}

	rownames(Theta)<-rownames(Y) 
	colnames(Theta)<-colnames(Y) 
	R <- Y - A - Theta

	addEffects <- list(mu = mu, Eeffect = Eeffect, Geffect = Geffect)
	coeff.list<-list(Y = meltName(Y,"Y"),
					 A = meltName(A,"A"), 
					 ThetaPlusR = meltName(Emat,"ThetaPlusR"), 
					 Theta = meltName(Theta,"Theta"), 
					 R = meltName(R,"R"))
	fitDF <- data.frame(Reduce(function(...) merge(..., by = c(E, G)), coeff.list), PCs)
	
	if (returnDataFrame){
		output <- list(DF = fitDF)
	} else {
		PClist <- sapply(colnames(PCs), function(x) NULL)
		for(i in 1:ncol(PCs)){
			PClist[[i]] <- matrix(PCs[,i], I, J)
		}
		output <- c(list(Y = Y, A = A, Theta = Theta, R = R), PClist)
	}

	forANOVA <- merge(dataReformatted[["DF"]], fitDF, all = TRUE, by = c(E, G))

	if (isUnRep) nPCANOVA <- KmaxPlusOne else nPCANOVA <- M
	modelStatement <- paste0(y, " ~ ", paste(c(E, G, colnames(PCs)[1:nPCANOVA]), collapse = " + "))

	if (!is.null(dataReformatted[["blockSig"]])){
		if (dataReformatted[["blockSig"]])  modelStatement <- paste0(modelStatement, " + ", block, ":", E)
	} 
	ANOVA <- as.data.frame(anova(lm(as.formula(modelStatement), data = forANOVA))) # use qr here instead of lm???
	names(ANOVA) <- c("Df", "SS", "MS", "Fval", "Pvalue")
	anovaDf <- degfGxE
	if (!dataReformatted[["isUnRep"]]) anovaDf <- c(anovaDf, tail(degfR, 1))
	ANOVA[grep("PC", rownames(ANOVA)), "Df"] <- anovaDf
	PCpvalue <- pvalue
	if (test == "bootstrap" & !isUnRep) PCpvalue <- c(PCpvalue, NA)

	ANOVA[nrow(ANOVA), "Df"] <- nrow(forANOVA) - sum(ANOVA[1:(nrow(ANOVA) - 1),"Df"]) - 1

	if(!is.null(errorMeanSqDfReps)) rownames(ANOVA)[nrow(ANOVA)] <- paste0("PC", M)

	
	ANOVA[grep("PC", rownames(ANOVA)), "Pvalue"] <- PCpvalue
	ANOVA[["MS"]] <- ANOVA[["SS"]] / ANOVA[["Df"]]
	ANOVA <- ANOVA[c("Df", "SS", "MS", "Pvalue")]
	if(!is.null(errorMeanSqDfReps)) ANOVA <- rbind(ANOVA, data.frame(Df = errorMeanSqDfReps[2], SS = prod(errorMeanSqDfReps[1:2]), MS = errorMeanSqDfReps[1], "Pvalue" = NA, row.names = "Residuals"))

	ANOVA[[" "]] <- stars(ANOVA[["Pvalue"]])
	
	ANOVA[["Pvalue"]] <- sprintf("%e", ANOVA[["Pvalue"]])
	ANOVA[ANOVA[["Pvalue"]] %in% "NA", "Pvalue"] <- NA

	ANOVA[as.numeric(ANOVA[["Pvalue"]]) %in% 0, "Pvalue"] <- paste0("< ", 1 / B)
	ANOVA[is.na(ANOVA[["Pvalue"]]), "Pvalue"] <- ""

	

	cat("Analysis of Variance Table\n", paste0("Response: ", y), "\n")
	print(ANOVA)
	cat("---\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05\n")

	cat("Number of significant multiplicative terms (tested sequentially): ", Kstar, "\n")
	options(contrasts = usrContr)
	return(c(model = model, addEffects, list(ANOVA = ANOVA), output, list(pvalue = pvalue, sigPC = Kstar, scores = scores, varcomp = Sigmasq,
	 		 degf = degf, svdE = Edecomp, rankTables = rankTables, nameRank = nameRank)))
}

