# Author: Nicholas Santantonio
# Institution: Cornell University, Plant Breeding and Genetics

# Version 1.8.5

# Please review the dimension selection method suggested by Forkman and Piepho (2014)
# the notation used within the script has been made to reflect thier terminology

# Forkman, J., & Piepho, H. P. (2014). Parametric bootstrap methods for testing multiplicative
# 		terms in GGE and AMMI models. Biometrics, 70(3), 639-647. 


bilinear <- function(x = NULL, G = NULL, E = NULL, y = NULL, block = NULL, model = "AMMI", errorMeanSqDfReps = NULL, f=0.5, test = "bootstrap", imputePC = "sig", alpha = 0.05, B = 100000, nCore = 1, Bonferroni = TRUE, returnDataFrame = TRUE, override3col = FALSE, ...){
# x = Y; f = 0.5; G = NULL; E = NULL; y = NULL; block = NULL; model = "AMMI"; test = "bootstrap"; errorMeanSqDfReps = NULL; alpha = 0.05; B = 10000; nCore = 2; Bonferroni = TRUE; returnDataFrame = TRUE; override3col = TRUE;
	list.of.packages <- c("reshape2","foreach", "doMC")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) {install.packages(new.packages)}
	
	packages.required<-c(TRUE, TRUE, nCore > 1)
	
	lapply(list.of.packages[packages.required], require, character.only = TRUE)

	if(nCore > 1) {registerDoMC(cores = nCore)}

	if(!.Platform$OS.type %in% "unix"){
		nCore <- 1
		warning("Multicore function only implemented for unix based systems (including OS X and Linux).
Continuing with 1 core... if analysis takes too long, press Ctrl+c to quit analysis.")
	}

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
	ANOVA <- as.data.frame(anova(lm(as.formula(modelStatement), data = forANOVA)))
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

AMMIplot <- function(bilinearObject, plots = "linear", color = c("darkgreen", "darkblue"),  PC = 1, f = 0.5, pdf.filename = NULL, ...){

	argumentChange <- function(defaultArgs, userArgs){
		userArgs <- list(...)
		defaultArgs[names(userArgs)] <- userArgs
		return(defaultArgs)
	}

	lineIntersect <- function(s1, i1, s2, i2){
		c( (i2 - i1) / (s1 - s2), (s1 * i2 - s2 * i1) / (s1 - s2))
	}

	Escores <- bilinearObject$scores$Escores
	Gscores <- bilinearObject$scores$Gscores

	Lambda <- bilinearObject$svdE$d
	mu <- bilinearObject$mu
	Geffect <- bilinearObject$Geffect
	Eeffect <- bilinearObject$Eeffect

	I <- length(Geffect)
	J <- length(Eeffect)
	
	M <- length(Lambda) - 1
	if(PC >= M) stop("Must specify 'PC' less than ", M)

	Kstar <- bilinearObject$sigPC
	model <- bilinearObject$model
	
	percExpl <- round(Lambda/sum(Lambda)*100,1)

	if(!all(plots %in% c("linear", "winner"))) stop("'plots' must be a character vector with arguments 'linear' and/or 'winner'")

	nplots <- length(plots) 
	if(nplots > 1) par(mfrow = c(1, nplots))

	Gintercept <- mu + Geffect
	Erange <- range(Escores[,paste0("PC", PC)])

	x <- seq(Erange[1], Erange[2], length.out = 1000)
	nominalLines <- sweep(x %*% t(Gscores[,paste0("PC", PC)]), 2, Gintercept, "+")
	winnerPos <- names(Geffect)[apply(nominalLines, 1, which.max)]
	winner <- unique(winnerPos)

	if("linear" %in% plots){
		labelx <- tapply(x, factor(winnerPos, levels = winner), mean)
		labely <- NULL

		for (i in 1:length(winner)){
			labely <- c(labely, Gintercept[winner[i]] + Gscores[winner[i],paste0("PC", PC)] * labelx[i])
		}
		labely <- labely + diff(range(nominalLines)) * 0.05

		xlimits <- range(Escores[,paste0("PC", PC)]) * 1.1 
		ylimits <- range(nominalLines - mu) * 1.1 + mu

		linecol <- rep("#000000", I)
		linecol[names(Geffect) %in% winner] <- color[1]

		linewd <- rep(1, I)
		linewd[names(Geffect) %in% winner] <- 2

		linPlotArgs <- argumentChange(list(main = "Linear AMMI plot", ylim = ylimits, xlim = xlimits, yaxs="i", pch = 17, 
							ylab = "Nominal", xlab = paste0("Environment PC", PC), cex = 1.5), list(...))

		if(!is.null(pdf.filename)) pdf(paste0("linear", pdf.filename))

		do.call(plot, c(list(x = Escores[,paste0("PC", PC)], y = rep(ylimits[1], dim(Escores)[1])), linPlotArgs))

		# linLineArgs <- list(col = linecol[i], lwd = linewd[i])

		for (i in 1:length(Gintercept)){
			# linLineDefaultArgsi <- list(col = linLineArgs$col[i], lwd = linLineArgs$lwd[i])
			# abline(Gintercept[i], Gscores[i,paste0("PC", PC)], linLineDefaultArgsi)
			abline(Gintercept[i], Gscores[i,paste0("PC", PC)], col = linecol[i], lwd = linewd[i])
		}
		text(labelx, labely, labels = winner, col = color[1])
		if(!is.null(pdf.filename)) dev.off()
	}

	if("winner" %in% plots){
		winnerInt <- list()
		for (i in 2:length(winner)){
			winnerInt[[i-1]] <- lineIntersect(Gscores[[winner[i-1],paste0("PC", PC)]], Gintercept[[winner[i-1]]], Gscores[[winner[i],paste0("PC", PC)]], Gintercept[[winner[i]]])
		}
		
		xlimits <- range(Eeffect) * 1.1 + mu
		ylimits <- range(Escores[, paste0("PC", PC)]) * 1.1
		
		if(!is.null(pdf.filename)) pdf(paste0("winners", pdf.filename))
		
		winPLotArgs <- argumentChange(list(xlim = xlimits, ylim = ylimits, main = "Winner AMMI plot", 
										   ylab = paste0("Environment PC", PC, "Score"), xlab = "Environmental mean", 
										   pch = 16, col = color[2]), list(...))
		do.call(plot, c(list(x = mu + Eeffect, y = Escores[,paste0("PC", PC)]), winPLotArgs))
		abline(h = sapply(winnerInt, function(x) x[1]))
		text(mu + Eeffect * 1.1, Escores[,paste0("PC", PC)] * 1.1, labels = rownames(Escores), col = color[2])
		
		winPointArgs <- argumentChange(list(pch = 16, col = color[1]), list(...))
		do.call(points, c(list(x = mu + Geffect [names(Geffect) %in% winner], y = Gscores[names(Geffect) %in% winner, paste0("PC", PC)]), winPointArgs))
		text(mu + Geffect [names(Geffect) %in% winner] * 1.1, Gscores[names(Geffect) %in% winner,paste0("PC", PC)] * 1.1, labels = rownames(Gscores)[names(Geffect) %in% winner], col = color[1])
		if(!is.null(pdf.filename)) dev.off()
	}
}


BBplot<-function(bilinearObject, f = 0.5, pdf.filename=NULL, nPC=NULL, plottitle="Biplot", biplot3D=FALSE, outer.box=TRUE, internal.axes=FALSE, scaled=TRUE, decorateGGE=FALSE, Gnames=TRUE){
	# many ideas for these plots from:  http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization

	Edecomp <- bilinearObject$svdE
	Lambda <- Edecomp$d
	Geffect <- bilinearObject$Geffect
	Eeffect <- bilinearObject$Eeffect
	M <- length(Lambda) - 1
	Kstar <- bilinearObject$sigPC
	model <- bilinearObject$model
	percExpl <- round(Lambda / sum(Lambda) * 100, 1)

	if(f < 0 | f > 1){ stop("'f' must be between 0 and 1, 
where f > 0.5 weights genotype scores more than environments,
and f < 0.5 weights environment scores more than genotype scores")}

	if(Kstar == 0){
		warning("No significant GxE detected! Biplots will only be produced if 'nPC' > 0.
Because the system appears to be completely additive, EXTREME CAUTION should be 
taken when interpreting any biplots produced from this data!")
	}

	if(is.null(nPC)){
		if(biplot3D) {maxPC <- 3} else {maxPC <- 2}
		if(model %in% "GGE" & !biplot3D){
			nPC <- 2
		} else {
			nPC<-min(Kstar,maxPC)
		}
	}

	if(nPC == 0){
		stop("Please specify 'nPC' as an integer > 0 to veiw biplots for NOT significant PCs.")
	}

	if(nPC > (M - 1)){
		stop("Data does not have enough dimensions for the number of specified PCs!")
	}

	if(Kstar < nPC){ warning("At least one of the plotted terms was found to be NOT significant in the previous analysis! 
Caution should be taken when evaluating the following biplot.")}

	if(biplot3D){
		list.of.packages <- c("rgl")
		new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
		if(length(new.packages)) {install.packages(new.packages)}
		
		packages.required<-c(biplot3D)
		lapply(list.of.packages[packages.required], require, character.only = TRUE)


		lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
		lim2 <- function(x){c(min(x), max(x)) * 1.1}

		rgl_add_axes <- function(x, y, z, axis.col = "grey10", xlab = "", ylab="", zlab="", show.plane = FALSE, show.bbox = TRUE, bbox.col = c("gray40","black")){ 
			xlim <- lim(x)
			ylim <- lim(y)
			zlim <- lim(z)
			rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
			rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
			rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
			  
			axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), c(0, 0, zlim[2]))
			rgl.points(axes, color = axis.col, size = 3)
			  
			rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col, adj = c(0.5, -0.8), size = 2)

			if(show.bbox){
				rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, emission=bbox.col[1], 
			          specular=bbox.col[1], shininess=5, xlen = 3, ylen = 3, zlen = 3) 
			  }
		}
	}

	if(biplot3D){
		if (nPC >= 3 & M >= 3){
			Gscores <- Edecomp$u[,(nPC-2):nPC] %*% diag(Lambda[(nPC-2):nPC]^(f))
			rownames(Gscores)<-names(Geffect)
			Escores <- Edecomp$v[,(nPC-2):nPC] %*% diag(Lambda[(nPC-2):nPC]^(1-f))
			rownames(Escores)<-names(Eeffect)
			axes.names<-paste(paste("PC",(nPC-2):nPC,sep=""), " ", percExpl[(nPC-2):nPC], "%", sep="")
			names(axes.names)<-c("labx","laby","labz")
		} else if(nPC < 3){
			Gscores <- cbind(Geffect,Edecomp$u[,1:2] %*% diag(Lambda[1:2]^(f)), deparse.level = 0)
			Escores <- cbind(Eeffect,Edecomp$v[,1:2] %*% diag(Lambda[1:2]^(1-f)), deparse.level = 0)
			axes.names<-c("Main Effect",paste(c("PC1","PC2"), " ", percExpl[1:2], "%", sep=""))
			names(axes.names)<-c("labx","laby","labz")
		} else {
			stop("'nPC' specified is larger than rank! Please specify less terms.")
		}


		if(!internal.axes){
			open3d()
    		plot3d(Gscores, cex=2, lwd=4, size = 7,box=outer.box,xlab="",ylab="",zlab="")
    		if(Gnames) text3d(Gscores * 1.1 , text=rownames(Gscores))
    		segments3d(x=c(t(cbind(0,Escores[,1]))),
       					y=c(t(cbind(0,Escores[,2]))),
       					z=c(t(cbind(0,Escores[,3]))), lwd=2, color="red")
    		text3d(Escores * 1.1 , text=rownames(Escores), color="red")
			if(scaled) {aspect3d(1,1,1)}
    		title3d(main=plottitle, xlab=axes.names["labx"],ylab=axes.names["laby"],zlab=axes.names["labz"])

			rgl.lines(lim(c(Gscores[,1],Escores[,1])), c(0, 0), c(0, 0), lwd=1, color = "gray15")
			rgl.lines(c(0, 0), lim(c(Gscores[,2],Escores[,2])), c(0, 0), lwd=1, color = "gray15")
			rgl.lines(c(0, 0), c(0, 0), lim(c(Gscores[,3],Escores[,3])), lwd=1, color = "gray15")

		} else {

			rgl.open()
			rgl.bg(color = "white")
			rgl.points(Gscores, cex=2, lwd=4, size = 7,color="black")
			text3d(Gscores*1.1, text=rownames(Gscores))
			segments3d(x=c(t(cbind(0,Escores[,1]))),
		       y=c(t(cbind(0,Escores[,2]))),
		       z=c(t(cbind(0,Escores[,3]))), lwd=2, color="red")
			text3d(Escores*1.1, text=rownames(Escores), color="red")
			if(scaled) {aspect3d(1,1,1)}
			
			rgl.lines(lim(c(Gscores[,1],Escores[,1])), c(0, 0), c(0, 0), lwd=1, color = "gray15")
			rgl.lines(c(0, 0), lim(c(Gscores[,2],Escores[,2])), c(0, 0), lwd=1, color = "gray15")
			rgl.lines(c(0, 0), c(0, 0), lim(c(Gscores[,3],Escores[,3])), lwd=1, color = "gray15")

			rgl_add_axes(x=c(Gscores[,1],Escores[,1]),
					y=c(Gscores[,2],Escores[,2]), 
					z=c(Gscores[,3],Escores[,3]),
					axis.col = "gray15",
					xlab=axes.names["labx"],
					ylab=axes.names["laby"],
					zlab=axes.names["labz"],
		            show.bbox = outer.box) 
			if(!outer.box){
				axis3d('x', pos=c( NA, 0, 0 ), col = "gray15")
				axis3d('y', pos=c( 0, NA, 0 ), col = "gray15")
				axis3d('z', pos=c( 0, 0, NA ), col = "gray15")
				title3d(main=plottitle)
    		} else {
    			title3d(main=plottitle, xlab=axes.names["labx"], ylab=axes.names["laby"], zlab=axes.names["labz"])
			}
		}

	} else {
		if(nPC == 1){
			Gscores <- cbind(Geffect,Lambda[1]^(f) * Edecomp$u[,1], deparse.level = 0)
			colnames(Gscores)<-c("Geffect","PC1")
			rownames(Gscores)<-names(Geffect)
			
			Escores <- cbind(Eeffect,Lambda[1]^(1-f) * Edecomp$v[,1], deparse.level = 0)
			colnames(Escores)<-c("Eeffect","PC1")
			rownames(Escores)<-names(Eeffect)
			
			axes.names<-c("Main Effect",paste("PC1", " ", percExpl[1], "%", sep=""))
			names(axes.names)<-c("labx","laby")
		} else {
			PCnames<-paste("PC",(nPC-1):nPC, sep="")
			
			Gscores <- Edecomp$u[,(nPC-1):nPC] %*% diag(Lambda[(nPC-1):nPC]^(f))
			colnames(Gscores)<-PCnames
			rownames(Gscores)<-names(Geffect)
			
			Escores <- Edecomp$v[,(nPC-1):nPC] %*% diag(Lambda[(nPC-1):nPC]^(1-f))
			colnames(Escores)<-PCnames
			rownames(Escores)<-names(Eeffect)
			
			axes.names<-paste(PCnames, " ", percExpl[(nPC-1):nPC], "%", sep="")
			names(axes.names)<-c("labx","laby")
		}
		if(!is.null(pdf.filename)){pdf(pdf.filename)}
		sect<-decorateBBplot(Gscores = Gscores, Escores = Escores, Eeffect = Eeffect, axes.names = axes.names, decorate = decorateGGE, scaledplot = scaled, Gnames=Gnames)
		if(!is.null(pdf.filename)){dev.off()}
		return(invisible(sect))
	}
}

svdSignFlipMC<-function(X, threshold=1e-15, nCore=1){

# method from:
# Bro, Rasmus, Evrim Acar, and Tamara G. Kolda. "Resolving the sign ambiguity in the 
# 		singular value decomposition." Journal of Chemometrics 22, no. 2 (2008): 135-140.
# Note: irrelavent when matrix is centered from both dimensions (?)

	MC <- nCore > 1
	
	list_of_packages <- c("foreach")
	if(MC) list_of_packages <- c(list_of_packages,"doMC")
	
	new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
	if(length(new_packages)) {install.packages(new_packages)}

	lapply(list_of_packages, require, character.only = TRUE)
	if(MC) registerDoMC(cores=nCore)

	UDV<-svd(X)
	K<-length(which(UDV$d > threshold))
	zeroeig <- length(UDV$d) - K

	signflip<-function(k, X, UDV){
		Dstar <- UDV$d
		Dstar[k] <- 0
		Y <- X - tcrossprod(UDV$u %*% diag(Dstar), UDV$v)
		uty <- crossprod(UDV$u[,k], Y)
		vty <- tcrossprod(UDV$v[,k], Y)

		sign_uty<-sign(uty)
		sign_vty<-sign(vty)

		l<- sum(sign_uty * uty^2)
		r<- sum(sign_vty * vty^2)

		if(l*r < 0){
			if(l < r){
				l <- l * -1
			} else {
				r <- r * -1
			}
		} 
		return(c(l,r))
	}

	if(MC){
		sf <- foreach(k = 1:K, .combine = rbind) %dopar% signflip(k,  X=X, UDV=UDV)
	} else {
		sf <- foreach(k = 1:K, .combine = rbind) %do% signflip(k,  X=X, UDV=UDV)
	}

	sl <- sf[,1]
	sr <- sf[,2]

	if(zeroeig > 0){
		sl <- c(sl, rep(1, zeroeig))
		sl <- c(sl, rep(1, zeroeig))
	}

	UDV$u <- sweep(UDV$u, 2, sign(sl), '*')
	UDV$v <- sweep(UDV$v, 2, sign(sr), '*')
	return(UDV)
}


decorateBBplot <- function(Gscores, Escores, Eeffect, axes.names, decorate = TRUE, scaledplot = TRUE, Gnames){

	lim2 <- function(x){c(min(x), max(x)) * 1.1}

	plot.labels<-function(mat,infl=0.01){
		plus<-apply(mat,2,sd) * infl 
		return(mat + sweep(sign(mat),2,plus,'*'))
	}

	Xprod2d <- function(a, b){
		det(rbind(a,b))
	}

	# essentially the same as xprod() by Grej https://gist.github.com/grej/3624583
	Xprod3d <- function(a, b){
		if (length(a) > 3 | length(a) > 3) { stop("vectors too long!")}

		if (length(a) < 3) { a <- c(a, rep(0, 3 - length(a)))}
		if (length(b) < 3) { b <- c(b, rep(0, 3 - length(b)))}
		
		ab<-rbind(a,b)
		x <- det(ab[,2:3])
		y <- -det(ab[,c(1,3)])
		z <- det(ab[,1:2])

		return(c(x,y,z))
	}

	XprodTest2d<-function(A,B,C){
		AxB<-Xprod2d(A,B)
		AxC<-Xprod2d(A,C)
		CxB<-Xprod2d(C,B)
		CxA<-Xprod2d(C,A)
		return(AxB * AxC >= 0 && CxB * CxA >=0)
	}

	if(scaledplot){
		limy <- limx <- lim2(c(Gscores,Escores))
	} else {
		limx <- lim2(c(Gscores[,1], Escores[,1]))
		limy <- lim2(c(Gscores[,2], Escores[,2]))
	}

	plot(0, type="n", xlab=axes.names["labx"], ylab=axes.names["laby"], xlim=limx, ylim=limy)
	axis(3, col="red")
	axis(4, col="red")
	
	chullxy <- chull(Gscores)
	if (decorate){

		outside <- rbind(Gscores[chullxy, ], Gscores[chullxy[1], ])

		slopes <- (diff(outside[,2])/diff(outside[,1]))
		nQuad <- length(slopes)

		intercepts<-Gscores[chullxy,2] - Gscores[chullxy,1] * slopes

		y<-NULL
		x<-NULL
		for (i in 1:length(intercepts)){
			A<-cbind(c(-slopes[i], 1/slopes[i]),1)
			b<-matrix(c(intercepts[i],0),ncol=1)
			xy<-solve(t(A) %*% A) %*% t(A) %*% b
			y[i]<-xy[2]
			x[i]<-xy[1]
		}

		y<-c(y,y[1])
		x<-c(x,x[1])
		Quads<-list()
		for(i in 1:nQuad){
			Quads[[i]]<-rbind(c(0,0),
							  c(x[i],y[i]),
							  outside[i+1,],
							  c(x[i+1],y[i+1]),
							  c(0,0)
				)
		}

		Evec<-as.list(data.frame(t(Escores)))

		sectors<-list()
		for (i in 2:length(x)){
			sectors[[i-1]]<-list(A=c(x[i-1], y[i-1]), C= c(x[i], y[i]))
		}
		names(sectors)<-paste("sector",1:length(sectors),sep="")


		ABC<-list()
		for (i in 1:length(Evec)){
			ABC[[i]]<-list()
			for (j in 1:length(sectors)){
				ABC[[i]][[j]]<-c(sectors[[j]], list(B=Evec[[i]]))
			}
			names(ABC[[i]])<-names(sectors)
		}
		names(ABC)<-names(Evec)

		inSector<-sapply(ABC,function(z){	which(sapply(z, function(zz) { XprodTest2d(A = zz$A, B = zz$B, C = zz$C) })) })
		megEnv<-sort(unique(inSector))

		megEnvEffect<-aggregate(Eeffect~inSector, FUN=mean, data=cbind(Eeffect,inSector))
		
		bluered<-colorRampPalette(c("#82AEC2","#BF4C3C"), alpha=TRUE)(nQuad)[rank(megEnvEffect["Eeffect"])]
		blueredlite30<-colorRampPalette(c("#82AEC24D","#BF4C3C4D"), alpha=TRUE)(nQuad)[rank(megEnvEffect["Eeffect"])]

		Ecols <- merge(data.frame(sector = inSector, E = gsub("\\.sector.*","",names(inSector)), ordered = 1:length(inSector)), data.frame(sector=1:length(bluered),color=bluered))
		Ecols <- Ecols[order(Ecols$ordered),]

		ii <- 0
		for (i in megEnv){
			ii <- ii + 1
			polygon(Quads[[i]],col=blueredlite30[ii], border=bluered[ii], lwd=2)
		}

		zero<-rep(0,length(megEnv) + 1)
		lineOrd<-unique(c(megEnv, megEnv + 1))
		lineOrd <- lineOrd[lineOrd %in% 1:nQuad]
		segments(zero, zero, x[lineOrd], y[lineOrd])

	}

	points(Gscores, pch=16)
	if(Gnames){
		text(plot.labels(Gscores,0.15),labels=rownames(Gscores))
	} else {
		text(plot.labels(Gscores,0.15)[chullxy,],labels=rownames(Gscores)[chullxy])
	}
	lines(x=c(t(cbind(0,Escores[,1]))), y=c(t(cbind(0,Escores[,2]))), lwd=2, col="red")

	text(plot.labels(Escores,0.15),labels=rownames(Escores), col="red")

	if(decorate) {
		return(invisible(list(Sector = Ecols)))
	}
}

