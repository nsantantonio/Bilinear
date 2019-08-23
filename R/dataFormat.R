#' dataFormat Function
#'
#' This function is made to format data within bilinear(). Not intended for use by user.
#'
#' @param ... arguments to dataFormat()
#' @return formatted data frame for use inside bilinear()
#' @keywords bilinear
#' @export
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
				cat("This program now handles missing cells in the genotype environment table! An EM algorithm will be used to impute missing cells\n")
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