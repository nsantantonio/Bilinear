#' BBplot function
#'
#' This function makes biplots from a bilinear() object, and is a wrapper for decorateBBplot.
#' 
#' @param bilinearObject object from output of a call to biliner()
#' @param f numeric. Scale parameter in (0, 1) for exponent on eigenvalues for weighting the genotype scores. Environment scores are weighted 1 - f. Default is 0.5 
#' @param nPC integer. Number of principal components to graph in the model. Default is 2 for a 2D biplot and 3 for a 3D biplot. If the number of significant PCs is less than the default and no argument is provided to 'nPC', the mean will be printed on the x axis, and the 1st (or 1st and 2nd for 3D biplot) PC will be printed on the remaining axis(es)
#' @param plottitle character. Title for biplot.
#' @param biplot3D logical. Should a 3dimensional 
#' @param outer.box logical. Shoudl a box be printed around the plot. default is TRUE.
#' @param internal.axes logical. For 3D biplots(?) default is FALSE 
#' @param scaled logical. default is TRUE
#' @param decorateGGE logical. should a GGE plot be decorate to dilineate megaenvironments to determine genotype winners?
#' @param Gnames character vector of genotype names to be included on plot near respective points
#' @details
#' put some details here
#'
#' @examples
#' 
#' data(soyMeanMat)
#' AMMIfit <- bilinear(x = soyMeanMat)
#' AMMIplot(AMMIfit)
#' AMMIplot(AMMIfit, "winner")
#' AMMIplot(AMMIfit, "winner", color = "hotpink")
#' AMMIplot(AMMIfit, c("linear", "winner"), color = c("hotpink", "darkorchid"))
#' 
#' data(ontario)
#' GGEfit <- bilinear(x = ontario, model = "GGE")
#' BBplot(GGEfit, decorateGGE = TRUE)
#'
#' @keywords AMMI, GGE, biplot
#' @importFrom utils installed.packages
#' @export
BBplot<-function(bilinearObject, f = 0.5, nPC=NULL, plottitle="Biplot", biplot3D=FALSE, outer.box=TRUE, internal.axes=FALSE, scaled=TRUE, decorateGGE=FALSE, Gnames=TRUE){
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
			nPC <- min(Kstar, maxPC)
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
		if(length(new.packages)) stop("Please install rgl to make 3D biplots: install.packages('rgl')")
		# if (!requireNamespace("rgl", quietly = TRUE)) stop("Please install rgl to make 3D biplots: install.packages('rgl')") # use this!
		
		packages.required <- c(biplot3D)
		lapply(list.of.packages[packages.required], require, character.only = TRUE)


		lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
		lim2 <- function(x){c(min(x), max(x)) * 1.1}

		rgl_add_axes <- function(x, y, z, axis.col = "grey10", xlab = "", ylab="", zlab="", show.plane = FALSE, show.bbox = TRUE, bbox.col = c("gray40","black")){ 
			xlim <- lim(x)
			ylim <- lim(y)
			zlim <- lim(z)
			rgl::rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
			rgl::rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
			rgl::rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
			  
			axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), c(0, 0, zlim[2]))
			rgl::rgl.points(axes, color = axis.col, size = 3)
			  
			rgl::rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col, adj = c(0.5, -0.8), size = 2)

			if(show.bbox){
				rgl::rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, emission=bbox.col[1], 
			          specular=bbox.col[1], shininess=5, xlen = 3, ylen = 3, zlen = 3) 
			  }
		}
	}

	if(biplot3D){

		# list_of_packages <- c("rgl")
		# new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

		# if(length(new_packages)){
		# 	cat("Additional Packages required to make 3D biplot:\n",list_of_packages, "\nDo you want to install them?\n")	
		# 	install.packages(new_packages)
		# }
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
    		plot3d(Gscores, cex = 2, lwd = 4, size = 7,box = outer.box, xlab = "", ylab = "", zlab =  "")
    		if(Gnames) text3d(Gscores * 1.1 , text = rownames(Gscores))
    		segments3d(x = c(t(cbind(0,Escores[, 1]))),
       					y = c(t(cbind(0,Escores[, 2]))),
       					z = c(t(cbind(0,Escores[, 3]))), lwd = 2, color = "red")
    		text3d(Escores * 1.1 , text = rownames(Escores), color = "red")
			if(scaled) {aspect3d(1, 1, 1)}
    		title3d(main = plottitle, xlab = axes.names["labx"], ylab = axes.names["laby"], zlab = axes.names["labz"])

			rgl.lines(lim(c(Gscores[, 1], Escores[, 1])), c(0, 0), c(0, 0), lwd = 1, color = "gray15")
			rgl.lines(c(0, 0), lim(c(Gscores[, 2], Escores[, 2])), c(0, 0), lwd = 1, color = "gray15")
			rgl.lines(c(0, 0), c(0, 0), lim(c(Gscores[, 3], Escores[, 3])), lwd = 1, color = "gray15")

		} else {

			rgl.open()
			rgl.bg(color = "white")
			rgl.points(Gscores, cex = 2, lwd = 4, size = 7, color = "black")
			text3d(Gscores * 1.1, text = rownames(Gscores))
			segments3d(x = c(t(cbind(0, Escores[, 1]))),
		       y = c(t(cbind(0, Escores[, 2]))),
		       z = c(t(cbind(0, Escores[, 3]))), lwd = 2, color = "red") 
			text3d(Escores * 1.1, text = rownames(Escores), color = "red")
			if(scaled) {aspect3d(1,1,1)}
			
			rgl.lines(lim(c(Gscores[, 1], Escores[, 1])), c(0, 0), c(0, 0), lwd = 1, color = "gray15")
			rgl.lines(c(0, 0), lim(c(Gscores[, 2], Escores[, 2])), c(0, 0), lwd = 1, color = "gray15")
			rgl.lines(c(0, 0), c(0, 0), lim(c(Gscores[, 3], Escores[, 3])), lwd = 1, color = "gray15")

			rgl_add_axes(x=c(Gscores[,1], Escores[,1]),
					y = c(Gscores[,2], Escores[,2]), 
					z = c(Gscores[,3], Escores[,3]),
					axis.col = "gray15",
					xlab = axes.names["labx"],
					ylab = axes.names["laby"],
					zlab = axes.names["labz"],
		            show.bbox = outer.box) 
			if(!outer.box){
				axis3d('x', pos = c( NA, 0, 0 ), col = "gray15")
				axis3d('y', pos = c( 0, NA, 0 ), col = "gray15")
				axis3d('z', pos = c( 0, 0, NA ), col = "gray15")
				title3d(main = plottitle)
    		} else {
    			title3d(main = plottitle, xlab = axes.names["labx"], ylab = axes.names["laby"], zlab = axes.names["labz"])
			}
		}

	} else {
		if(nPC == 1){
			Gscores <- cbind(Geffect, Lambda[1]^(f) * Edecomp$u[,1], deparse.level = 0)
			colnames(Gscores) <- c("Geffect","PC1")
			rownames(Gscores) <- names(Geffect)
			
			Escores <- cbind(Eeffect, Lambda[1]^(1-f) * Edecomp$v[,1], deparse.level = 0)
			colnames(Escores) <- c("Eeffect","PC1")
			rownames(Escores) <- names(Eeffect)
			
			axes.names <- c("Main Effect", paste("PC1", " ", percExpl[1], "%", sep=""))
			names(axes.names) <- c("labx", "laby")
		} else {
			PCnames<-paste("PC",(nPC-1):nPC, sep="")
			
			Gscores <- Edecomp$u[,(nPC-1):nPC] %*% diag(Lambda[(nPC-1):nPC]^(f))
			colnames(Gscores) <- PCnames
			rownames(Gscores) <- names(Geffect)
			
			Escores <- Edecomp$v[,(nPC-1):nPC] %*% diag(Lambda[(nPC-1):nPC]^(1-f))
			colnames(Escores) <- PCnames
			rownames(Escores) <- names(Eeffect)
			
			axes.names <- paste(PCnames, " ", percExpl[(nPC-1):nPC], "%", sep="")
			names(axes.names) <- c("labx","laby")
		}
		sect <- decorateBBplot(Gscores = Gscores, Escores = Escores, Eeffect = Eeffect, axes.names = axes.names, decorate = decorateGGE, scaledPlot = scaled, Gnames=Gnames)
		return(invisible(sect))
	}
}
