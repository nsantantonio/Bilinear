# Author: Nicholas Santantonio
# Institution: Cornell University, Plant Breeding and Genetics

source("~/Google Drive/Bilinear/Bilinear v1.6.R") # this will have to be changed to the directory where Bilinear v1.6.R lives on your computer!

# Please review the dimension selection method suggested by Forkman and Piepho (2014)
# the notation used within the script has been made to reflect thier terminology

# Forkman, J., & Piepho, H. P. (2014). Parametric bootstrap methods for testing multiplicative
# 		terms in GGE and AMMI models. Biometrics, 70(3), 639-647. 

# input New York soybean data used by Forkman and Piepho (2014)
Environment <- c("A77", "A86", "C87", "C88", "G85", "G88", "I85", "N87", "R81","V79") 
Genotype <- c("CHIP", "CORS", "EVAN", "HODG", "S200", "WELL", "WILK")

Y <- matrix(c(
2333.25, 2815.00, 2326.25, 3097.00, 2812.75,
3509.75, 1587.50, 2529.75, 2350.25, 1278.00,
3106.50, 3644.00, 1715.50, 2822.25, 3192.75,
4209.00, 1798.25, 2703.00, 2967.50, 1661.00, 
2725.25, 3164.00, 3273.00, 3923.50, 3257.75, 
4906.25, 1735.50, 2543.50, 2037.50, 1111.25, 
2741.00, 3407.25, 2573.50, 3164.50, 3422.75, 
4751.75, 1839.25, 2830.75, 2800.50, 1499.25, 
2843.00, 3188.25, 1237.75, 2482.00, 3207.25, 
4501.50, 1374.75, 2583.50, 2617.00, 1963.25, 
2745.25, 3589.75, 1660.75, 2901.50, 2925.25, 
3689.75, 1533.75, 2538.25, 2758.00, 1711.75, 
2470.50, 2605.75, 3141.25, 3612.50, 2960.50, 
4608.25, 1606.50, 2284.50, 1386.50, 578.00), 
nrow = 10, ncol = 7,
dimnames = list(Environment, Genotype))
print(Y)

# typically data is structured such that Genotypes live in rows and Environments live in columns, 
# but this is arguably debatable, and really depends on how you are veiwing the problem.
# this function REQUIRES (kinda...) that genotypes be on rows and environments on columns for GGE/SREG
# or GREG models, however it shouldnt matter for AMMI model fit, but will affect the biplots 
# developed from these fits in that genotypes will instead be shown as vectors and environments 
# as points. The data must also contain at least 3 genotypes measured in at least 3 environments.
Y<-t(Y) 
print(Y)  

#################
# The functions #
#################

# bilinear() function to fit bilinear models and test for significant multiplicative terms
# bilinear(x=NULL,G=NULL,E=NULL,y=NULL, model="AMMI", alpha=0.05, B=100000,
#		 nCore=1, Bonferroni=TRUE, returnDataFrame=TRUE, pdf.filename=NULL)

# BBplot() function to produce biplots from a bilinear model fit object
# BBplot(bilinearObject, pdf.filename=NULL, nPC=NULL, plottitle="Biplot", 
# 	   biplot3D=FALSE, outer.box=FALSE, internal.axes=FALSE, scaled=TRUE)


###########################################################
# Additive Main effects Multiplicative Interaction (AMMI) #
###########################################################

# fit AMMI model and test for significant multiplicative terms,  
# where y = mu + E_i + G_j + sum(u_k * d_k * v_k) + e_ijkl for k = 1,...,K multiplicative terms 
# Note: environment and genotype indices are left out for simplicity from here on
AMMIfit <- bilinear(x = Y)  # note this may take several minutes, as it is sampling 100,000 bootstraps with 1 cpu

# Check p-values relative to those obtained by Forkman and Piepho (2014).  While they will likely not match exactly, they should be VERY close.

##############
### Output ###
##############
# bilinear returns a list of relevant objects to the model fit

# main effects
AMMIfit[c("mu", "E_effect", "G_effect")] 

# data.frame in long format of model coefficients from model using the notation of Forkman and Piepho (2014)
# y = mu + E + G + sum(u_k * d_k * v_k) + e for k = 1,...,K significant multiplicative terms
# where A = y - sum(u_k * d_k * v_k) + e, such that A is the sum of additive effects, and in the AMMI case, A = mu + G + E
# Theta = sum(u_k * d_k * v_k) for k = 1,...,K significant multiplicative terms and is the GxE deviations
# R = e, i.e. the residuals
# and ThetaPlusR = Theta + R. This is called E by Forkman and Piepho (2014), 
# and is the residual matrix from the additive model, but is renamed here to avoid confusion with the 'E' for environment.
# the model can therefore be rewritten in matrix format as Y = A + Theta + R
# If returnDataFrame = FALSE, then the coefficients will be returned as I x J matrices for I genotypes in J environments
AMMIfit["df"] 

# P-values and the number of 'significant' terms determined by a sequential hypothesis test
# of multiplicative terms until the null hypothesis is not rejected.
# by default, a multiple test correction is obtained by dividing alpha by K for the kth term.
# ie. the 1st term is deemed significant if the parametric bootstrap pvalue < alpha for the first term,
# pvalue < alpha/2 for the second term (given that the first hypothesis was rejected), pvalue < alpha/3 
# (given the first two hypotheses were rejected),  ... , and so on. 
# However, Forkman and Peipho seem to suggest that this is unnecessary due to the sequential testing method
# and the multiple test correction can be overridden by setting Bonferroni = FALSE in the bilinear() function call
AMMIfit[c("Pval","sigPC")] # Note that the first two terms are considered significant

# Estimates of GxE and error (across environments) variance components computed as shown by Forkman and Peipho (2014)
# If a balanced set of phenotypic observations within environment are provided, the within error variance is also
# calculated and returned with the Sigma vector.
AMMIfit["Sigma"]

# singular value decomposition of the ThetaPlusR matrix (called E by Forkman and Piepho (2014))
# Note that due to floating point algebra, the last eigenvalue is not zero, but is essentially so.
# this lack of rank is the result of subtracting the additive effects (A) from Y, and restricting
# the rows and columns to sum to zero. The rank of ThetaPlusR (i.e. E), is therefore min(I-1,J-1) for the AMMI,
# min(I,J-1) for the GGE/SREG, and min(I-1,J) for the GREG.
AMMIfit["svd_E"]


#################################
# Some model fit considerations #
#################################

# The number of bootstraps ('B' argument) can be reduced to say 10,000 to reduce run time.  
# For more accurate P-values, higher bootstrap numbers are required. 
# the authors suggest the number of bootstraps be on the order of 1 x 10^5
# however, 1 x 10^4 or even 1 x 10^3 might be used if computation time becomes a significant problem 
# this is unlikely to be an issue for typical numbers Genotypes and Environments, because it is limited by the smaller of the two.
system.time(AMMIfit <- bilinear(x = Y, B = 10000))

# use multiple cores ('nCore' argument) to reduce run time by approximately 1/min(nCore, J-1, I-1), where the GE matrix is of dimension IxJ 
# this is limited by the number of CPUs your computer has (most computers have at least 2 these days),
# and may not work outside of the terminal environment (e.g. not tested in Rstudio or R GUI).
# additionally, the multicore functionallity is limited to macOS, Linux and most other unix based operating systems, 
# and will not work with Microsoft Windows (because 'dmMC' package doesnt play well with Windows).
system.time(AMMIfit <- bilinear(x = Y, model = "AMMI", B = 10000, nCore = 2))

# Smoke em if you got em!!! 
# system.time(AMMIfit<-bilinear(x=Y, model="AMMI", B=100000)) # slow
# system.time(AMMIfit<-bilinear(x=Y, model="AMMI", B=100000, nCore=5)) # much faster!!!!

##################################
# Biplots of bilinear model fits #
##################################

# default 2D biplot of significant terms from AMMI fit
BBplot(AMMIfit) 

# write 2D biplot to PDF (note: 3D plots cannot be saved as PDF)
BBplot(AMMIfit, pdf.filename = "testBiplot.pdf")
# BBplot(AMMIfit, pdf.filename = "AMMIexample2PC.pdf") # Ignore....
# BBplot(AMMIfit, scaled = TRUE, pdf.filename = "AMMIexample2PCscaled.pdf") # Ignore...

# If there are many genotypes, can use Gnames = FALSE to only display the outer hull, or "ring", of genotypes
# this has yet to be implemented in the 3D version, and setting Gnames = FALSE will only display points
BBplot(AMMIfit, Gnames = FALSE) 

# Typically the parameter f is 0.5, such that genotype and environment scores are weighted equally, 
# however this can be changed such that f > 0.5 weights genotype scores more than environments,
# and f < 0.5 weights environment scores more than genotype scores
BBplot(AMMIfit, f = 0.6) 
BBplot(AMMIfit, f = 0.4)

# 3D biplot of significant terms and main effects (i.e. genotype and environment means)
BBplot(AMMIfit, biplot3D = TRUE)
BBplot(AMMIfit, biplot3D = TRUE, Gnames = FALSE)


# remove outer box
BBplot(AMMIfit, biplot3D = TRUE, outer.box = FALSE)

# unscaled biplot, might be useful to interpret scale, but visually unappealing.  Other plots have proper scaling on each dimension,
# but the axes are just streched to make a square box
BBplot(AMMIfit, biplot3D = TRUE, scaled = FALSE)

# Slightly different 3D biplot of significant terms and main effects 
BBplot(AMMIfit, biplot3D = TRUE, internal.axes = TRUE)
BBplot(AMMIfit, biplot3D = TRUE, outer.box = FALSE, internal.axes = TRUE)


# 2D biplot of main effect and first term
# the 'nPC' argument will override the significance of multiplicative terms, 
# plotting the nPC-1 x nPC terms for biplot3D=FALSE and the nPC-2 x nPC-1 x nPC terms for biplot3D=TRUE 
# otherwise, the appropriate plot should be produced automatically.
BBplot(AMMIfit, nPC = 1)

#Oops! That was hard to see because the main effect is MUCH larger than the interaction effect, lets not scale it
BBplot(AMMIfit, nPC = 1, scaled = FALSE)
# BBplot(AMMIfit, nPC =1, scaled = FALSE, pdf.filename = "AMMIexample1PC.pdf") # Ignore...

# 3D biplot of second and third terms, regardless of 3rd term significance 
BBplot(AMMIfit, nPC = 3)

# 3D biplot of first three terms, regardless of 3rd term significance 
BBplot(AMMIfit, nPC = 3, biplot3D = TRUE)


################
# Other Models #
################


# GGE or Sites Regression (SREG) for 'who wins where'. However, the AMMI will typically 
# show this information in a 2D biplot if nPC=1 or in a 3D biplot if nPC=2. This is because 
# in SREG, the first term is typically directly related to the additive effect of Genotypes.

# GGE fit, where y = mu + E + sum(u_k * d_k * v_k) + e for k = 1,...,K multiplicative terms
# 'model' can also be specified as "SREG"
GGEfit<-bilinear(x = Y, model = "GGE", B = 10000, nCore = 2)
BBplot(GGEfit)

# decorateGGE can only be used with the GGE model type, and draws "mega-Environments" see Yan et al 2007
BBplot(GGEfit, decorateGGE = TRUE)
# BBplot(GGEfit, decorateGGE=TRUE, pdf.filename="GGEexample2PC.pdf")
BBplot(GGEfit, biplot3D = TRUE)

# Genotypes regression model (GREG). Used to diagnose environments. Similarly to the SREG, the AMMI
# will typically show this in a 2D biplot if nPC=1 or in a 3D biplot if nPC=2
# GREG fit, where y = mu + G + sum(u_k * d_k * v_k) + e for k = 1,...,K multiplicative terms
GREGfit<-bilinear(x = Y, model = "GREG", B = 10000, nCore = 2)
BBplot(GREGfit)
BBplot(GREGfit, biplot3D = TRUE) 

#################################
# Additional ways to input data #
#################################

# the bilinear function is relatively flexible to several different data types in 
# long or wide formats.  Typically the data is submitted as a dataframe (wide or long)
# or matrix of Genotype within Environment means. Data MUST BE BALANCED, with all 
# genotypes observed in all environments to fit this model! 

# Additionally if data is submitted as a dataframe of raw values the function will 
# calculate means within each environment and run the analysis on the resulting 
# means table. A within environment error term is also estimated and returned in 
# the output. This is completely equivalent to submitting a table of means if the
# data is perfectly balanced (even for an RCBD). However, results may be ambiguous
# if the data is unbalanced.

# type bilinear() in the terminal for a more detailed description of data formats

#reformat data
require(reshape2)
dfY <- melt(Y)

# as long as dataframe or matrix column names match c("G", "E", "y") 
# exactly (order doesnt matter), the function will assume that the data is in long format
names(dfY) <- c("G","E","y") # name dfY c("G", "E", "y") 
print(head(dfY, 20))
bilinear(x = dfY, model = "AMMI", B = 10000, nCore = 2)

matY <- as.matrix(dfY)
print(matY[1:20,])
bilinear(x = matY, model = "AMMI", B = 10000, nCore = 2) 

# if the columns do not match as above, the names of the variables must be
#  passed to the 'G', 'E', and 'y' arguments as character strings
names(dfY) <- c("geno", "env", "trait") #change names dfY
print(head(dfY, 20))
bilinear(x = dfY, G = "geno", E = "env", y = "trait", model = "AMMI", B = 10000, nCore = 2)

# data can be also supplied directly as vectors to the 'G', 'E' and 'y' arguments
Gvec <- dfY$geno
Evec <- dfY$env
yvec <- dfY$trait
bilinear(G = Gvec, E = Evec, y = yvec, model = "AMMI", alpha = 0.05, B = 10000, nCore = 2)





#######################################################################
# Additional Example using Onterio Wheat Data in                      #
# 	Yan, W., Hunt, L. A., Sheng, Q., & Szlavnics, Z. (2000).          #
#		Cultivar evaluation and mega-environment investigation based  #
#		on the GGE biplot. Crop Science, 40(3), 597-605.              #
#######################################################################

Onterio<-matrix(c(
4.46,4.15,2.85,3.08,5.94,4.45,4.35,4.04,2.67,
4.42,4.77,2.91,3.51,5.70,5.15,4.96,4.39,2.94,
4.67,4.58,3.10,3.46,6.07,5.03,4.73,3.90,2.62,
4.73,4.75,3.38,3.90,6.22,5.34,4.23,4.89,3.45,
4.39,4.60,3.51,3.85,5.77,5.42,5.15,4.10,2.83,
5.18,4.48,2.99,3.77,6.58,5.05,3.99,4.27,2.78,
3.38,4.18,2.74,3.16,5.34,4.27,4.16,4.06,2.03,
4.85,4.66,4.43,3.95,5.54,5.83,4.17,5.06,3.57,
5.04,4.74,3.51,3.44,5.96,4.86,4.98,4.51,2.86,
5.20,4.66,3.60,3.76,5.94,5.35,3.90,4.45,3.30,
4.29,4.53,2.76,3.42,6.14,5.25,4.86,4.14,3.15,
3.15,3.04,2.39,2.35,4.23,4.26,3.38,4.07,2.10,
4.10,3.88,2.30,3.72,4.56,5.15,2.60,4.96,2.89,
3.34,3.85,2.42,2.78,4.63,5.09,3.28,3.92,2.56,
4.38,4.70,3.66,3.59,6.19,5.14,3.93,4.21,2.93,
4.94,4.70,2.95,3.90,6.06,5.33,4.30,4.30,3.03,
3.79,4.97,3.38,3.35,4.77,5.30,4.32,4.86,3.38,
4.24,4.65,3.61,3.91,6.64,4.83,5.01,4.36,3.11),
 ncol=9, nrow=18, byrow=TRUE, dimnames=list(paste("G",1:18,sep=""), paste("E",1:9,sep="")))

GGEfit<-bilinear(Onterio, nCore= 2 , model = "GGE", B = 10000)
BBplot(GGEfit, decorateGGE = TRUE)



