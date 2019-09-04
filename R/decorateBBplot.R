#' decorateBBplot Function
#'
#' This function is called by BBplot for drawing GGE biplot with environment winner designations. This function is not intended to be used by the user directly
#'
#' @param Gscores matrix of genotype scores
#' @param Escores matrix of environment scores
#' @param Eeffect numeric vector environment effects(?)
#' @param axes.names character vector of axis names
#' @param decorate logical, indicates if the environment winenr panels should be colored
#' @param scaledPlot logical, indicates if genotypes scores and escores shoudl be scaled to one another (how does f effect this?)
#' @param Gnames names of genotypes to be printed, default is outer edge
#' @return list sectors
#' @details 
#' 
#' Uses chull() to compute convex hull of biplot points, draws lines to each of them sucessively, and then draws lines from the origin perpendicular to each segment.
#' Colors are added to sectors containing genotypes if decorate = TRUE
#'
#' @examples
#' 
#' data(ontario)
#' GGEfit <- bilinear(x = ontario, model = "GGE")
#' BBplot(GGEfit, decorateGGE = TRUE)
#'
#' @keywords GGE 
#' @importFrom grDevices chull colorRampPalette 
#' @importFrom graphics axis lines plot points polygon segments text
#' @importFrom stats aggregate
#' @export
decorateBBplot <- function(Gscores, Escores, Eeffect, axes.names, decorate = TRUE, scaledPlot = TRUE, Gnames){

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

	if(scaledPlot){
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
			A <- cbind(c(-slopes[i], 1 / slopes[i]), 1)
			b <- matrix(c(intercepts[i], 0), ncol = 1)
			# xy <- solve(crossprod(A) %*% crossprod(A, b))
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
							  outside[i + 1,],
							  c(x[i + 1], y[i + 1]),
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
