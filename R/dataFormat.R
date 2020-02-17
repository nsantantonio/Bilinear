#' dataFormat Function
#'
#' This function is made to format data within bilinear(). Not intended for use by user.
#'
#' @param x matrix or data.frame containing numeric values of cell means with genotypes on rows and environments on columns (e.g. ith genotype within jth environment). 
#'			Alternatively a data.frame can be provided in long format, with factor and value variable names passed to 'G', 'E' and 'y' respectively
#' @param G character. Name of genotype variable in data.frame 'x' if 'x' is a data.frame in long format. Optionally, character or factor vector of genotype names (optional). If NULL (default) a matrix of cell means must be provided to 'x' argument
#' @param E character. Name of environment variable in data.frame 'x' if 'x' is a data.frame in long format. Optionally, character or factor vector of environment names (optional). If NULL (default) a matrix of cell means must be provided to 'x' argument
#' @param y character. Name of phenotype variable in data.frame 'x' if 'x' is a data.frame in long format. Optionally, numeric vector of phenotype values (optional). If NULL (default) a matrix of cell means must be provided to 'x' argument
#' @param block character. Optional, for RCBD only. Name of block variable in data.frame 'x' if 'x' is a data.frame in long format and study is an RCBD. Optionally, character or factor vector of block names. 
#' @param alpha pvalue cutoff threshold for significance. Default is 0.05.
#' @param anyNullGEy logical. Are any of 'G', 'E', 'y' arguments NULL?
#' @param override3col logical. Overrides the 3 column numeric warning when a matrix of 3 coulmns is provided. Most users will not need to turn this off, unless the number of environments is 3 and a table of cell means are provided. Default is FALSE.
#' @return formatted data frame for use inside bilinear()
#' @keywords bilinear
#' @importFrom stats anova as.formula dummy.coef lm
#' @export
dataFormat <- function(x, G, E, y, block, alpha, anyNullGEy, override3col){
	meltName <- function(x, G, E, vName) {
		melted <- data.frame(rep(rownames(x), ncol(x)), rep(colnames(x), each = nrow(x)), c(x))
		names(melted) = c(G, E, vName)
		melted
	}
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
			if (!override3col) is3env <- readline("You have provided a dataframe or matrix with 3 columns, and no variable names to 'G', 'E' and/or 'y'.
						Is the provided dataframe or matrix contain 3 environments with genotypes on the rows?  (y/n) ")
			if (override3col | is3env %in% c("y","yes","YES","Yes","Y")) {
				x <- as.matrix(x)
			} else {
				stop(warnmessage)
			}
		}

		if (all(c("G","E","y") %in% colnames(x))){
			DF <- as.data.frame(x[,c("G", "E", "y")])
			G <- "G"; E <- "E";	y <- "y"
			if (!class(DF[[y]]) %in% "numeric") DF[[y]] <- as.numeric(DF[[y]])
			if (isRCBD) {
				DF$block <- x[,"block"]
				block <- "block"
			}
		} else  if (is.data.frame(x)){
			DF <- x
		} else if (is.matrix(x)){
			if (is.numeric(x) & allNullGEy) {
				if(is.null(rownames(x))) rownames(x) <- paste0("G", 1:nrow(x))
				if(is.null(colnames(x))) colnames(x) <- paste0("E", 1:ncol(x))
				DF <- meltName(x, G = "G", E = "E", vName = "y")
				return(list(Y = x, DF = DF, isUnRep = TRUE, sigmasq = NULL, repPerG = 1)) #, anyMissCells = anyMissCells))
			} else if (is.numeric(x) & !allNullGEy){
				DF <- x
			} else if(!is.numeric(x)) {
				stop("The matrix supplied is of class: ", class(c(x)), ". \nPlease provide a matrix with numeric values (i.e. is.numeric(x) should be TRUE)")
			} else {
				stop("Data format failed. please submit an issue with a reproducable example at https://github.com/nsantantonio/Bilinear/issues")
			}
		} else {
			stop(warnmessage)
		}			
	}

	if (any(sapply(DF, class) == "factor")) DF <- droplevels(DF)
	
	if (is.integer(DF[[E]]) | suppressWarnings(sum(is.na(as.numeric(as.character(DF[[E]]))))) == 0) DF[[E]] <- as.character(DF[[E]])
	if (is.integer(DF[[G]]) | suppressWarnings(sum(is.na(as.numeric(as.character(DF[[G]]))))) == 0) DF[[G]] <- as.character(DF[[G]])

	Elvls <- if (is.factor(DF[[E]])) levels(DF[[E]]) else unique(DF[[E]][!is.na(DF[[E]])])
	Glvls <- if (is.factor(DF[[G]])) levels(DF[[G]]) else unique(DF[[G]][!is.na(DF[[G]])])


	if (!is.factor(DF[[E]])) DF[[E]] <- factor(DF[[E]], levels = Elvls)
	if (!is.factor(DF[[G]])) DF[[G]] <- factor(DF[[G]], levels = Glvls)

	repGE <- as.matrix(table(DF[[G]][!is.na(DF[[y]])], DF[[E]][!is.na(DF[[y]])]))
	repGE[Glvls, Elvls]
	repGE[, Elvls]

	nReps <- unique(repGE)
	isUnRep <- all(nReps[nReps != 0] == 1)
	
	I <- nrow(repGE)
	J <- ncol(repGE)

	if (I < 3 | J < 3) {stop(warnlackrank)}
	
	if (length(nReps) > 1){ if (length(nReps[nReps != 0]) > 1) cat(warnunbal) } 

	if (isUnRep){
		fit <- lm(as.formula(paste0(y," ~ ", E, " + ", G)), data = DF)
		allCombos <- expand.grid(Elvls, Glvls)
		names(allCombos) <- c(E, G)
		DFwNA <- merge(allCombos, DF, by = c(E, G), all = TRUE)
		Y <- matrix(DFwNA[[y]], I, J, dimnames = list(unique(DFwNA[[G]]), unique(DFwNA[[E]])))
		Y <- Y[Glvls, Elvls]
	} else {
		if (isRCBD){	
			fit <- lm(as.formula(paste0(y," ~ ", E, " + ", block , ":", E, " + ", G, " + ", G, ":", E)), data = DF)
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
		intercept <- coefs[["(Intercept)"]]
		whichInt <- which(c(paste0(G, ":", E), paste0(E, ":", G)) %in% names(coefs))

		Emat <- matrix(rep(coefs[[E]], each = I), nrow = I, ncol = J)
		Gmat <- matrix(rep(coefs[[G]], J), nrow = I, ncol = J)
		GEcoefs <- coefs[names(coefs) %in% c(paste0(G, ":", E), paste0(E, ":", G))][[1]]

		allCombos <- expand.grid(E = Elvls, G = Glvls)		
		colnames(allCombos) <- if (whichInt == 2) c(E, G) else c(G, E) 
		if (!all(names(GEcoefs) == apply(allCombos, 1, paste, collapse = ":"))) stop("G:E terms are in the wrong order! Please post an issue at https://github.com/nsantantonio/Bilinear/issues")
		GEmat <- matrix(GEcoefs, nrow = I, ncol = J, byrow = whichInt == 2)
		
		Y <- intercept + Gmat + Emat + GEmat
		rownames(Y) <- Glvls
		colnames(Y) <- Elvls
		Y[repGE == 0] <- NA
	}
	if (isUnRep) sigmasq <- NULL else sigmasq <- summary(fit)$sigma^2
	errorDf <- anova(fit)["Residuals", "Df"]
	returnList <- list(Y = Y, DF = DF, isUnRep = isUnRep)

	if (!isUnRep){
		sigmaList <- c(sigmasq = sigmasq, errorDf = errorDf)
		returnList <- c(returnList, list(sigmasq = sigmaList))
	}
	if (length(nReps) > 1) repPerG <- repGE else repPerG <- nReps
	if (isRCBD) returnList <- c(returnList, blockSig = blockSig)
	returnList <- c(returnList, list(repPerG = repPerG))
	return(returnList)
}