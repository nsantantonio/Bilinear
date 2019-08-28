# Author: Nicholas Santantonio
# Institution: Cornell University, Plant Breeding and Genetics

library(Bilinear) # load bilinear package


# Please review the dimension selection method suggested by Forkman and Piepho (2014)
# the notation used within the script has been made to reflect thier terminology

# Forkman, J., & Piepho, H. P. (2014). Parametric bootstrap methods for testing multiplicative
# 		terms in GGE and AMMI models. Biometrics, 70(3), 639-647. 

# load example dataset from Zobel, R. W., Wright, M. J., & Gauch, H. G. (1988). Statistical analysis of a yield trial. Agronomy Journal, 80(3), 388-393.
data(soy)
print(soyMeanMat)

# typically data is structured such that Genotypes live in rows and Environments live in columns, 
# but this is arguably debatable, and really depends on how you are veiwing the problem.
# this function REQUIRES (kinda...) that genotypes be on rows and environments on columns for GGE/SREG
# or GREG models, however it shouldnt matter for AMMI model fit, but will affect the biplots 
# developed from these fits in that genotypes will instead be shown as vectors and environments 
# as points. The data must also contain at least 3 genotypes measured in at least 3 environments.

#################
# The functions #
#################

# documentation
# ?bilinear
# ?BBplot
# ?AMMIplot


###########################################################
# Additive Main effects Multiplicative Interaction (AMMI) #
###########################################################

# fit AMMI model and test for significant multiplicative terms,  
# where y = mu + E_i + G_j + sum(u_k * d_k * v_k) + e_ijkl for k = 1,...,K multiplicative terms 
# Note: environment and genotype indices are left out for simplicity from here on
AMMIfit <- bilinear(x = soyMeanMat)  # note this could take several minutes, as it is sampling 100,000 bootstraps with 1 cpu
# Check p-values relative to those obtained by Forkman and Piepho (2014).  While they will likely not match exactly, they should be VERY close.

##############
### Output ###
##############
# bilinear returns a list of relevant objects to the model fit

# main effects
AMMIfit[c("mu", "Eeffect", "Geffect")] 

# data.frame in long format of model coefficients from model using the notation of Forkman and Piepho (2014)
# y = mu + E + G + sum(u_k * d_k * v_k) + e for k = 1,...,K significant multiplicative terms
# where A = y - sum(u_k * d_k * v_k) + e, such that A is the sum of additive effects, and in the AMMI case, A = mu + G + E
# Theta = sum(u_k * d_k * v_k) for k = 1,...,K significant multiplicative terms and is the GxE deviations
# R = e, i.e. the residuals
# and ThetaPlusR = Theta + R. This is called E by Forkman and Piepho (2014), 
# and is the residual matrix from the additive model, but is renamed here to avoid confusion with the 'E' for environment.
# the model can therefore be rewritten in matrix format as Y = A + Theta + R
# If returnDataFrame = FALSE, then the coefficients will be returned as I x J matrices for I genotypes in J environments
AMMIfit["DF"] 

# P-values and the number of 'significant' terms determined by a sequential hypothesis test
# of multiplicative terms until the null hypothesis is not rejected.
# by default, a multiple test correction is obtained by dividing alpha by K for the kth term.
# ie. the 1st term is deemed significant if the parametric bootstrap pvalue < alpha for the first term,
# pvalue < alpha/2 for the second term (given that the first hypothesis was rejected), pvalue < alpha/3 
# (given the first two hypotheses were rejected),  ... , and so on. 
# However, Forkman and Peipho seem to suggest that this is unnecessary due to the sequential testing method
# and the multiple test correction can be overridden by setting Bonferroni = FALSE in the bilinear() function call
AMMIfit[c("pvalue", "sigPC")] # Note that the first two terms are considered significant

# Estimates of GxE and error (across environments) variance components computed as shown by Forkman and Peipho (2014)
# If a balanced set of phenotypic observations within environment are provided, the within error variance is also
# calculated and returned with the Sigma vector.
AMMIfit["varcomp"]

# singular value decomposition of the ThetaPlusR matrix (called E by Forkman and Piepho (2014))
# Note that due to floating point algebra, the last eigenvalue is not zero, but is essentially so.
# this lack of rank is the result of subtracting the additive effects (A) from Y, and restricting
# the rows and columns to sum to zero. The rank of ThetaPlusR (i.e. E), is therefore min(I-1,J-1) for the AMMI,
# min(I,J-1) for the GGE/SREG, and min(I-1,J) for the GREG.
AMMIfit["svdE"]

#################################
# Some model fit considerations #
#################################

# The number of bootstraps ('B' argument) can be reduced to say 1,000 to reduce run time.  
# For more accurate P-values, higher bootstrap numbers are required. 
# the authors suggest the number of bootstraps be on the order of 1 x 10^5
# however, 1 x 10^4 or even 1 x 10^3 might be used if computation time becomes a significant problem 
# this is unlikely to be an issue for typical numbers Genotypes and Environments, because it is limited by the smaller of the two.
system.time(AMMIfit <- bilinear(x = soyMeanMat, B = 1000))

# use multiple cores ('nCore' argument) to reduce run time by approximately 1/min(nCore, J-1, I-1), where the GE matrix is of dimension IxJ 
# this is limited by the number of CPUs your computer has (most computers have at least 2 these days),
# and may not work outside of the terminal environment (e.g. not tested in Rstudio or R GUI).
# Additionally, the multicore functionalty needs a parallel backend for 'foreach', such as 'doMC'.
library(doMC)
registerDoMC(2) # number of cores
system.time(AMMIfit <- bilinear(x = soyMeanMat, model = "AMMI", B = 10000, nCore = 2))
# note that here, parallelization actually slows down the process because the job is small, and the overhead actually increases computation time
# with larger problems and more bootstraps, parallelization should drastically decrease computation time. 
# system.time(AMMIfit<-bilinear(x=soyMeanMat, model="AMMI", B=100000, nCore=5)) # much faster!!!!


############################
# F tests for significance #
############################

# since the data is simply cell means and we have no way to calculate a within environment error term, 
# 'Ftest' tests using a sequential F-test known to be too liberal, see below for F_R test for replicated data
bilinear(soyMeanMat, test = "Ftest")[c("pvalue", "sigPC")]

###################################
# Raw Replicates and Missing Data # 
# Using same soybean data from    #
# Gauch 1988					  #
###################################

soyboot <- bilinear(soy, B = 10000)
# Here we will use the Fr test from Piepho (1995) to test for significant PCs
soyF <- bilinear(soy, test = "Ftest")

# We can also use the residual error variance estimate to refit with 
# genotype within environment means .
# This is useful if you have means and a single error variance

soyfit <- lm(as.formula(y ~ E + G + G:E + block:E), data = soy)
errorMeanSq <- summary(soyfit)$sigma^2
errorDf <- anova(soyfit)["Residuals", "Df"]
reps <- unique(table(soy$G, soy$E))

print(c(errorMeanSq, errorDf, reps))
soyF2 <- bilinear(soyMeanMat, errorMeanSqDfReps = c(errorMeanSq, errorDf, reps), test = "Ftest")


# small example with no sig GxE
pea <- matrix(c(3505, 3779, 3951, 3826, 3997,
 4081, 4207, 4278, 4065, 4274,
 4238, 4371, 4378, 4315, 4552,
 2844, 3203, 3345, 3224, 3193),
 ncol = 5, nrow = 4, byrow = TRUE, dimnames  = list(paste0("E", 1:4), paste0("G", 1:5)))
AMMIpea <- bilinear(x = t(pea)[,1:3], override3col = TRUE, B = 10000)

###########################
# linear and winner plots #
###########################

# For Details about these plots see Figures 2 and 3 of Gauch, H., & Zobel, R. W. (1997). Identifying mega-environments and targeting genotypes. Crop Science, 37(2), 311-326.

AMMIplot(AMMIfit)
AMMIplot(AMMIfit, "winner")
AMMIplot(AMMIfit, "winner", color = "hotpink")
AMMIplot(AMMIfit, c("linear", "winner"), color = c("hotpink", "darkorchid"))

dev.off()
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
GGEfit <- bilinear(x = soyMeanMat, model = "GGE", B = 10000)
BBplot(GGEfit)

# decorateGGE can only be used with the GGE model type, and draws "mega-Environments" see Yan et al 2007
BBplot(GGEfit, decorateGGE = TRUE)
# BBplot(GGEfit, decorateGGE=TRUE, pdf.filename="GGEexample2PC.pdf")
BBplot(GGEfit, biplot3D = TRUE)

# Genotypes regression model (GREG). Used to diagnose environments. Similarly to the SREG, the AMMI
# will typically show this in a 2D biplot if nPC=1 or in a 3D biplot if nPC=2
# GREG fit, where y = mu + G + sum(u_k * d_k * v_k) + e for k = 1,...,K multiplicative terms
GREGfit <- bilinear(x = soyMeanMat, model = "GREG", B = 10000)
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

# as long as dataframe or matrix column names match c("G", "E", "y") 
# exactly (order doesnt matter), the function will assume that the data is in long format
print(head(soyMeanDf))
bilinear(x = soyMeanDf, model = "AMMI", B = 10000)


matY <- as.matrix(dfY)
print(matY[1:20,])
bilinear(x = matY, model = "AMMI", B = 10000) 

# if the columns do not match as above, the names of the variables must be
#  passed to the 'G', 'E', and 'y' arguments as character strings

soyNewNames <- soy
names(soyNewNames) <- c("geno", "env", "block", "trait") #change names dfY
print(head(dfY, 20))
bilinear(x = dfY, G = "geno", E = "env", y = "trait", block = "block", model = "AMMI", B = 10000)

# data can be also supplied directly as vectors to the 'G', 'E' and 'y' arguments
Gvec <- dfY$geno
Evec <- dfY$env
yvec <- dfY$trait
bilinear(G = Gvec, E = Evec, y = yvec, model = "AMMI", alpha = 0.05, B = 10000)



# with missing cells; DOES NOT WORK YET!!!!
# bilinear(soyMiss10, test = "Ftest")


#######################################################################
# Additional Example using Onterio Wheat Data in                      #
# 	Yan, W., Hunt, L. A., Sheng, Q., & Szlavnics, Z. (2000).          #
#		Cultivar evaluation and mega-environment investigation based  #
#		on the GGE biplot. Crop Science, 40(3), 597-605.              #
#######################################################################
data(ontario)


GGEfit<-bilinear(ontario, model = "GGE", B = 10000)
BBplot(GGEfit, decorateGGE = TRUE)
