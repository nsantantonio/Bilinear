#' AMMIplot function
#'
#' This function will make linear and winner plots from a bilinear() object like Figure 2 and Figure 3 of Gauch & Zobel (1997).
#'
#' @param bilinearObject object from output of a call to biliner()
#' @param plots character vector of maximum length 2 containing the names of the plots to be produced. Possible arguments are "linear" and "winner". If both are specified in a single call, both plots will be plotted to the same device
#' @param color character vector of length 2 containing colors for the plots. The first specifies the genotype color and the second the environment color. If only one color is specified, only genotypes will be colored
#' @param PC integer. Principal component to plot. The default is 1
#' @param f numeric. Scale parameter in (0, 1) for exponent on eigenvalues for weighting the genotype scores. Environment scores are weighted 1 - f. Default is 0.5 
#' @param pdf.filename name of pdf file to output. 
#' @details
#' Many arguments can be passed on to plot() through (...)
#'
#' @examples
#' 
#' data(soy)
#' AMMIfit <- bilinear(x = soyShort)
#' AMMIplot(AMMIfit)
#' AMMIplot(AMMIfit, "winner")
#' AMMIplot(AMMIfit, "winner", color = "hotpink")
#' AMMIplot(AMMIfit, c("linear", "winner"), color = c("hotpink", "darkorchid"))
#' 
#' @keywords AMMI
#' @export

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


