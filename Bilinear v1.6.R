# Author: Nicholas Santantonio
# Institution: Cornell University, Plant Breeding and Genetics

# Version 1.6

# Please review the dimension selection method suggested by Forkman and Piepho (2014)
# the notation used within the script has been made to reflect thier terminology

# Forkman, J., & Piepho, H. P. (2014). Parametric bootstrap methods for testing multiplicative
# 		terms in GGE and AMMI models. Biometrics, 70(3), 639-647. 

bilinear<-function(x = NULL, G = NULL, E = NULL, y = NULL, model = "AMMI", alpha = 0.05, B = 100000, nCore = 1, Bonferroni = TRUE, returnDataFrame = TRUE, signflip = FALSE, MultiCorePackage = "doMC"){

	list.of.packages <- c("reshape2","foreach", MultiCorePackage)
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) {install.packages(new.packages)}
	
	packages.required<-c(TRUE, TRUE, nCore > 1)
	
	lapply(list.of.packages[packages.required], require, character.only = TRUE)

	if(nCore>1 & MultiCorePackage %in% "doMC") {registerDoMC(cores=nCore)}

	if(!.Platform$OS.type %in% "unix"){
		nCore <- 1
		warning("Multicore function not implemented for non - unix based systems (including OS X and Linux).
Continuing with 1 core... if analysis takes too long, press ctrl+c to quit analysis.")
	}

	meltName<-function(x,v.name){
		x<-melt(x)
		names(x)<-c("G","E",v.name)
		return(x)
	}

	bootstrap<-function(K, D, Dtilde, Lambda, B){
		Kplus1 <- K + 1
		Lambda_k <- Lambda[(Kplus1) : length(Lambda)]	
		T <- c(Lambda[Kplus1]^2 / (t(Lambda_k) %*% Lambda_k))

		T_b <- matrix(NA,nrow = B, ncol = 1)
		b <- 1
		while(b <= B){
			Ehatb<-matrix(rnorm((D - K) * (Dtilde - K)),  nrow = (D - K), ncol = (Dtilde - K))  # switched row and column dimensions... is this right?
			Lambda_b <- svd(Ehatb)$d
			T_b[b,1] <- Lambda_b[1]^2 / crossprod(Lambda_b) 
			b <- b+1
		}

		Pval<-sum(T_b > T)/B
		return(Pval)
	}

	extractPC<-function(k,UDV){
		d <- UDV$d
		d[!(1:length(d) %in% k)] <- 0
		PC <- c(UDV$u %*% diag(d) %*% t(UDV$v))
		return(PC)
	}



	warnmessage<-c("Please provide a dataframe or matrix for x!
If in wide format, put genotypes as rows and environments
as columns. If in long format, provide the names of the dataframe 
variables that correspond to genotypes, environments and the trait value 
as character strings to the 'G', 'E' and 'y' arguments, respectively.

Alternatively, you may provide vectors of genotype names, 
environment names and the trait value directly to the 'G', 
'E' and 'y' arguments, respectively.\n")

	warnlackrank<-"At least 3 Genotypes and 3 Environments are necessary for these bilinear models!\n"

	warnunbal<-"All genotypes must be observed at least once in every environment! 
Unbalanced designs are not recommended, but means within sites can be
used if replications differ. Care should be taken when drawing inference
from results of unbalanced designs."

	anynullGEy <- is.null(G) | is.null(E) | is.null(y)
	isallnum <- all(sapply(x,class) %in% c("numeric", "integer")) 
	isallch <- !any(is.numeric(G) | is.numeric(E) | is.numeric(y)) & is.null(x)

	
	if(is.null(x)){
		if (anynullGEy | isallch) {
			stop(warnmessage)
		} else {
			Glev <- length(unique(G))
			Elev <- length(unique(E))
			if(Glev < 3 | Elev < 3) {stop(warnlackrank)}
			Ncells <- Glev * Elev

			df<-data.frame(G,E,y)

			if(Ncells < nrow(df)){
				sigma.withinEnv <- summary(lm(as.formula(y ~ G + E)))$sigma
				df <- aggregate(as.formula(y ~ G + E), FUN = mean, data = df)
			} else if(Ncells > nrow(df)){
				stop(warnunbal)
			} 
			Y<-acast(df, G~E, value.var="y")
		}
	} else if (is.matrix(x)){
		if (all(colnames(x) %in% c("G","E","y"))){
				Y<-acast(as.data.frame(x), as.formula(paste("G ~ E")), value.var="y")
				class(Y)<-"numeric"
			} else{
				Y<-x
			}
		if(any(dim(Y) < 3)) {stop(warnlackrank)}
	} else if (ncol(x)==3){
		if (all(colnames(x) %in% c("G","E","y"))){
				Y<-acast(x, G~E, value.var="y")
		} else if(anynullGEy){
				is.3env<-readline("You have provided a dataframe with 3 columns, and no variable names to 'G', 'E' and/or 'y'.
						Is the provided dataframe a matrix of 3 environments with genotypes on the rows?  (y/n) ")
				if(is.3env %in% c("y","yes","YES","Yes","Y")) {
					Y<-as.matrix(df)
				} else {
					stop(warnmessage)
				} 
		} else {
			Glev <- length(unique(x[[G]]))  
			Elev <- length(unique(x[[E]]))
			if(Glev < 3 | Elev < 3) {stop(warnlackrank)}
			Ncells <- Glev * Elev

			if(Ncells < nrow(x)){
				sigma.withinEnv<-summary(lm(as.formula(paste(y, "~", G,"+", E)), data = x))$sigma
				df<-aggregate(as.formula(paste(y, "~", G,"+", E)), FUN = mean, data = x)
			} else if(Ncells > nrow(x)){
				stop(warnunbal)
			}  else {
				df<-x
			}
			Y<-acast(df, as.formula(paste(G,"~",E)), value.var=y)
		}
	} else if(isallnum){
		if(any(dim(Y) < 3)) {stop(warnlackrank)}
		Y<-as.matrix(x)
	} else {
		stop(warnmessage)
	}

	cat("\nDisplaying first (up to) 10 rows and columns of two way table of genotype within environment means:\n\n")
	print(Y[1:min(nrow(Y), 10), 1:min(ncol(Y), 10)])
	cat("\n")

	I<-nrow(Y) 
	J<-ncol(Y) 

	cat("\nData consists of ",I," Genotypes evaluated in ",J," Environments\n\n")

	if(model == "AMMI"){
		LLt_I <- diag(I) - (1 / I) * matrix(1, I, I) 
		LLt_J <- diag(J) - (1 / J) * matrix(1, J, J) 
		M<-min(I - 1, J - 1)
		D <- I - 1
		Dtilde <- J - 1
		E <- LLt_I %*% Y %*% LLt_J
	} else if (model %in% c("GGE","SREG")){
		LLt_I <- diag(I) - (1 / I) * matrix(1, I, I) 
		M <- min(I-1, J)
		D <- I - 1 
		Dtilde <- J
		E <- LLt_I %*% Y
	} else if (model %in% c("EGE","GREG")){
		LLt_J <- diag(J) - (1 / J) * matrix(1, J, J) 
		M <- min(I, J-1)
		D <- I - 1
		Dtilde <- J
		E <- Y %*% LLt_J
	 	# Need to check this to make sure correct, as  Forkman and Piepho dont cover GREG, 99% sure its right...
	} else {
		stop("'model' argument must be either 'AMMI', 'GGE', 'SREG', 'EGE',or 'GREG'!
Note: 'GGE' and 'SREG' are equivalent, as are 'EGE' and 'GREG'.")
	}
	rownames(E)<-rownames(Y) 
	colnames(E)<-colnames(Y) 
	nu <- D*Dtilde
	
	Kmax <- M - 2
	A <- Y - E
	
	if(signflip){
		Edecomp <- svdSignFlipMC(E,nCore = nCore)
	} else {
		Edecomp <- svd(E)
	}
	
	Lambda <- Edecomp$d

	if(nCore>1){
		PCs <- foreach(k=1:length(Edecomp$d), .combine=cbind) %dopar% extractPC(k, UDV = Edecomp)
	} else {
		PCs <- foreach(k=1:length(Edecomp$d), .combine=cbind) %do% extractPC(k, UDV = Edecomp)
	}
	colnames(PCs) <- paste0("PC",1:ncol(PCs))

	# AHHHHH Still need to figure this out!!!!!!!!!!!!!! SVD sign ambiguity
	# ThetaPlusR <- meltName(E,"ThetaPlusR")
	# PCdf <-  cbind(ThetaPlusR, PCs)
	# PCbeta <- lm(as.formula(paste0("ThetaPlusR ~ ", paste(colnames(PCs[,1:(ncol(PCs) - 1)]), collapse = "+"))), data=PCdf)$coefficients
	# PCsign <- sign(PCbeta[!names(PCbeta) %in% "(Intercept)"])

	if(nCore>1){
		Pval<-foreach(K= 0:Kmax, .combine = c) %dopar% bootstrap(K, D = D, Dtilde = Dtilde, Lambda = Lambda, B = B)
	} else {
		Pval<-foreach(K= 0:Kmax, .combine = c) %do% bootstrap(K, D = D, Dtilde = Dtilde, Lambda = Lambda, B = B)
	}

	if(Bonferroni){
		alpha <- alpha / 1:(Kmax + 1)	
	}
	names(Pval)<- paste("term", 1:length(Pval),sep="")
	sig <- Pval < alpha

	Kstar<-min(which(!sig)) - 1

	cat("P-values of multiplicative terms: \n"); print(Pval); cat("\n")
	cat("Number of significant multiplicative terms (tested sequentially): ", Kstar, "\n")

	if( Kstar > 1){
		sigmasqGxEhat<- 1/nu * crossprod(Lambda[1:Kstar])
		Theta<-Edecomp$u[,1:Kstar] %*% tcrossprod(diag(Lambda[1:Kstar]), Edecomp$v[,1:Kstar])	
	} else if( Kstar == 1){
		sigmasqGxEhat<- 1/nu * Lambda[1]^2
		Theta<-Lambda[1] * tcrossprod(Edecomp$u[,1], Edecomp$v[,1])
	} else {
		sigmasqGxEhat<- 0
		Theta<-matrix(0,nrow=nrow(Y),ncol=ncol(Y))
	}

	sigmasqhat <- 1/nu * crossprod(Lambda[(Kstar+1) : length(Lambda)])
	Sigma<-sqrt(c(sigmaGxE=sigmasqGxEhat,sigma=sigmasqhat))

	if(exists("sigma.withinEnv")) {Sigma["sigma.withinEnv"] <- sigma.withinEnv}

	rownames(Theta)<-rownames(Y) 
	colnames(Theta)<-colnames(Y) 
	R <- Y - A - Theta

	mu<- c("(Intercept)" = mean(Y))
		E_effect = colMeans(Y) - mu
		G_effect = rowMeans(Y) - mu

	if(returnDataFrame){
		add.effects<-list(mu = mu, E_effect = E_effect, G_effect = G_effect)

		coeff.list<-list(Y = meltName(Y,"Y"),
							A = meltName(A,"A"), 
							ThetaPlusR = meltName(E,"ThetaPlusR"), 
							Theta = meltName(Theta,"Theta"), 
							R = meltName(R,"R"))
		output <- c(model=model, add.effects, list(df = data.frame(Reduce(function(...) merge(..., by=c("E","G")), coeff.list), PCs), Pval=Pval, sigPC=Kstar, Sigma=Sigma, svd_E=Edecomp))
	} else {
		output <- c(model=model, add.effects, list(Y=Y, A=A, Theta=Theta, R=R, Pval=Pval, sigPC=Kstar, Sigma=Sigma), svd_E=Edecomp, list(PCs = PCs))
	}

	return(output)
}


BBplot<-function(bilinearObject, f=0.5, pdf.filename=NULL, nPC=NULL, plottitle="Biplot", biplot3D=FALSE, outer.box=TRUE, internal.axes=FALSE, scaled=TRUE, decorateGGE=FALSE, Gnames=TRUE){
	# many ideas for these plots from:  http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization

	Edecomp <- bilinearObject$svd_E
	Lambda <- Edecomp$d
	G_effect <- bilinearObject$G_effect
	E_effect <- bilinearObject$E_effect
	M <- length(Lambda) - 1
	Kstar <- bilinearObject$sigPC
	model <- bilinearObject$model
	perc.expl<-round(Lambda/sum(Lambda)*100,1)

	if(f < 0 | f > 1){ stop("'f' must be between 0 and 1, 
where f > 0.5 weights genotype scores more than environments,
and f < 0.5 weights environment scores more than genotype scores")}

	if(Kstar==0){
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

	if(nPC==0){
		stop("Please specify 'nPC' as an integer > 0 to veiw biplots for NOT significant PCs.")
	}

	if(nPC>(M-1)){
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
			rownames(Gscores)<-names(G_effect)
			Escores <- Edecomp$v[,(nPC-2):nPC] %*% diag(Lambda[(nPC-2):nPC]^(1-f))
			rownames(Escores)<-names(E_effect)
			axes.names<-paste(paste("PC",(nPC-2):nPC,sep=""), " ", perc.expl[(nPC-2):nPC], "%", sep="")
			names(axes.names)<-c("labx","laby","labz")
		} else if(nPC < 3){
			Gscores <- cbind(G_effect,Edecomp$u[,1:2] %*% diag(Lambda[1:2]^(f)), deparse.level = 0)
			Escores <- cbind(E_effect,Edecomp$v[,1:2] %*% diag(Lambda[1:2]^(1-f)), deparse.level = 0)
			axes.names<-c("Main Effect",paste(c("PC1","PC2"), " ", perc.expl[1:2], "%", sep=""))
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
    			title3d(main=plottitle, xlab=axes.names["labx"],ylab=axes.names["laby"],zlab=axes.names["labz"])
			}
		}

	} else {
		if(nPC == 1){
			Gscores <- cbind(G_effect,Lambda[1]^(f) * Edecomp$u[,1], deparse.level = 0)
			colnames(Gscores)<-c("G_effect","PC1")
			rownames(Gscores)<-names(G_effect)
			
			Escores <- cbind(E_effect,Lambda[1]^(1-f) * Edecomp$v[,1], deparse.level = 0)
			colnames(Escores)<-c("E_effect","PC1")
			rownames(Escores)<-names(E_effect)
			
			axes.names<-c("Main Effect",paste("PC1", " ", perc.expl[1], "%", sep=""))
			names(axes.names)<-c("labx","laby")
		} else {
			PCnames<-paste("PC",(nPC-1):nPC, sep="")
			
			Gscores <- Edecomp$u[,(nPC-1):nPC] %*% diag(Lambda[(nPC-1):nPC]^(f))
			colnames(Gscores)<-PCnames
			rownames(Gscores)<-names(G_effect)
			
			Escores <- Edecomp$v[,(nPC-1):nPC] %*% diag(Lambda[(nPC-1):nPC]^(1-f))
			colnames(Escores)<-PCnames
			rownames(Escores)<-names(E_effect)
			
			axes.names<-paste(PCnames, " ", perc.expl[(nPC-1):nPC], "%", sep="")
			names(axes.names)<-c("labx","laby")
		}
		if(!is.null(pdf.filename)){pdf(pdf.filename)}
		sect<-decorateBBplot(Gscores = Gscores, Escores = Escores, E_effect = E_effect, axes.names = axes.names, decorate = decorateGGE, scaledplot = scaled, Gnames=Gnames)
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


decorateBBplot <- function(Gscores, Escores, E_effect, axes.names, decorate = TRUE, scaledplot = TRUE, Gnames){

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

		megEnvEffect<-aggregate(E_effect~inSector, FUN=mean, data=cbind(E_effect,inSector))
		
		bluered<-colorRampPalette(c("#82AEC2","#BF4C3C"), alpha=TRUE)(nQuad)[rank(megEnvEffect["E_effect"])]
		blueredlite30<-colorRampPalette(c("#82AEC24D","#BF4C3C4D"), alpha=TRUE)(nQuad)[rank(megEnvEffect["E_effect"])]

		Ecols<-merge(data.frame(sector=inSector, E=gsub("\\.sector.*","",names(inSector)), ordered=1:length(inSector)), data.frame(sector=1:length(bluered),color=bluered))
		Ecols<-Ecols[order(Ecols$ordered),]

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

