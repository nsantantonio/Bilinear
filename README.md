# ```Bilinear``` package for fitting bilinear models (AMMI, GGE, etc.)

# v0.2 release

# Package formatting

<!-- (August 21, 2017, eclipse edition)  This is the first formal R package for ```bilinear```, version 0.1!  To install, first install ```devtools``` by entering an R environment and typing.  -->

To install, first install ```devtools``` by entering an R environment and typing.  

```install.packages("devtools")```

Load the devtools library with

```library(devtools)```

Then install bilinear by typing

```install_github("nsantantonio/Bilinear")```

You can now load the library as you would any other R package

```library(Bilinear)```

I highly reccomend that you do not download these files directly! Follow the installation instructions above.

# Microsoft R open

I have had issues using devtools with Microsoft R open <https://mran.microsoft.com/open>. I have supplied a source tarball, ```Bilinear_0.2.0.tar.gz``` for those who use MRAN. You can also just clone the repo and build on your own machine. 

```git clone https://github.com/nsantantonio/Bilinear```
```R CMD build Bilinear```
```R CMD INSTALL Bilinear```

# Missing cells now allowed!
Missing genotype/environment combinations are now allowed, and are imputed with an expectation maximization algorithm within ```bilinear()```. There is also functionality to allow the user to determine the accuracy of imputation, using the ```em()``` function (although this function is not required to be called by the user under typical circumstances).

# ```Bilinear```

Fit bilinear models using AMMI (Additive Main effects Multiplicative Interaction) or GGE/SREG (Genotype and Genotype by Environment /Sites REGression) and make 2D and 3D biplots.  Also tests for the number of significant dimensions of GxE using a parameteric bootstrap test method suggested by Forkman and Piepho (2014). An additional test for significant multiplicative terms, the F<sub>R</sub> test from Piepho (1995), has been also been implemented in version 0.1.

type ```?bilinear``` for some examples.  ```Bilinear``` version 0.1 contains the relevant functions, namely ```bilinear()``` to fit model and run statistical tests for signficant dimensions, ```BBplot()``` to plot 2D and 3D biplots, and ```AMMIplot()``` to produce linear and winner plots for AMMI models such as those in Gauch & Zobel (1997).

### Plant Breeders Disclaimer
It is this author's opinion that GGE type models are more difficult to interpret than AMMI type models, as the genotypic main effect will be present in one or more of the remaining dimensions of the residual matrix, but not necessarily the first (largest) dimension (e.g. when Var(GxE) > Var(G)).  Therefore, I recommend using AMMI type models, as both genotype and environment main effects are removed from the residuals before GxE effects are assessed.  However, GGE type models are provided here as a freely available resource for those who want to use them. The ```decorateGGE``` option of BBplot will draw the mega-environment delineations suggested by Yan et al. (2000).

## Proceed with Caution
This software is in development and bugs may (almost certainly) be present. I have tried to make the package flexible to different data inputs, but this leads to difficulties in testing all the possible ways one might try to input data. See known bugs below.

## Known issues that have been fixed

* Unbalanced data is now allowed! Missing genotype/environment cells are imputed using an expectation maximization algorithm. 

* ANOVA table for replicated data is now correct. The SS and MS were previously incorrectly calculated for the ANOVA table when the data were replicated, but the model fit, effects, p-values, and other results were unaffected. 

## Known issues to be addressed

* genotype and environment names CANNOT be integers, and must contain at least one non-numeric character (it should also work if they are of class ```character```, but I have yet to test this). This will be fixed in the very near future when I get a few spare minutes to look into it. 

* The plotting functions will be updated to allow more flexibility to the user to add certain line names, change colors etc.
<!-- * the dataframe needs to be sorted by environment and then genotype in order to print the proper mean square values for the PCs. This does not effect the test, but it does make an erroneous Mean Square values for the PCs.  -->

<!-- * There is an inconsistancy in the behavior of the winner plot of ```AMMIplot()``` with some data. I have not had the time to investigate this and do not know when I will get to it. -->


<!-- * Need balanced data.  At this time, the program requires balanced data across genotypes and environments (i.e. all genotypes observed in all environments).  If just one or two cells are missing, you could impute the genotype effect + environment effect (i.e. no GxE) for that cell. The program will run with unequal replication within each location, but each genotype must be observed at least once in each environment, and unequal replication could result in erroneous estimates (the program should print a warning if there is unequal replication).  Eventually an EM algorithm might be implemented to account for unbalanced data as suggested by Gauch and Zorbel (1990). This is planned to be implemented soon. For unbalanced data, where each genotype is observed at least once in each environment, you could cetainly use a mixed model treating genotypes as random and estimate within environment blups.   -->

<!-- ### Eventually R package to CRAN ?:
When I find time to produce a slightly more polished, flexible and tested program, it may be submitted to CRAN as an R package. 
 -->
### Soy data
The ```soy.Rdata``` included here is from Zobel, Wright and Gauch (1988) and is subset from the ```gauch.soy``` dataset included in the [```agridat``` package](https://github.com/kwstat/agridat). This was included to show the equivalence between using within environment means and raw data. 
### Ontario data
The ```ontario.Rdata``` included here is the Onterio wheat data containing 18 genotypes evaluated in 9 environments from from Yan et al. (2000) transcribed directly from the Yan et al. (2007) paper. 


## References 
##### Imputation using expectation maximization
- Gauch, H. G., & Zobel, R. W. (1990). Imputing missing yield trial data. Theoretical and Applied Genetics, 79(6), 753-761.

##### Bootstrap test
- Forkman, J., & Piepho, H. P. (2014). Parametric bootstrap methods for testing multiplicative terms in GGE and AMMI models. Biometrics, 70(3), 639-647. 

##### F test
- Piepho, H.P., 1995. Robustness of statistical tests for multiplicative terms in the additive main effects and multiplicative interaction model for cultivar trials. Theoretical and Applied Genetics, 90(3-4), pp.438-443.

##### SVD sign ambiguity (might not be relevant for centered matrices...)
- Bro, Rasmus, Evrim Acar, and Tamara G. Kolda. "Resolving the sign ambiguity in the singular value decomposition." Journal of Chemometrics 22, no. 2 (2008): 135-140.

#### selected background

##### AMMI
- Zobel, R. W., Wright, M. J., & Gauch, H. G. (1988). Statistical analysis of a yield trial. Agronomy Journal, 80(3), 388-393.
- Gauch, H., & Zobel, R. W. (1997). Identifying mega-environments and targeting genotypes. Crop Science, 37(2), 311-326.
- Gauch, H. G. (2013). A simple protocol for AMMI analysis of yield trials. Crop Science, 53(5), 1860-1869.

##### GGE
- Yan, W., Hunt, L. A., Sheng, Q., & Szlavnics, Z. (2000). Cultivar evaluation and mega-environment investigation based on the GGE biplot. Crop Science, 40(3), 597-605.

##### EM for unbalanced data
- Gauch Jr, H. G., & Zobel, R. W. (1990). Imputing missing yield trial data. Theoretical and Applied Genetics, 79(6), 753-761.

##### Discussion of bilinear models and biplot interpretation. I highly recommend reading them in the order shown here:

- Gauch, H. G. (2006). Statistical analysis of yield trials by AMMI and GGE. Crop science, 46(4), 1488-1500.
- Yan, W., Kang, M. S., Ma, B., Woods, S., & Cornelius, P. L. (2007). GGE biplot vs. AMMI analysis of genotype-by-environment data. Crop science, 47(2), 643-653.
- Gauch, H. G., Piepho, H. P., & Annicchiarico, P. (2008). Statistical analysis of yield trials by AMMI and GGE: Further considerations. Crop Science, 48(3), 866-889.
- Yang, R. C., Crossa, J., Cornelius, P. L., & Burgueño, J. (2009). Biplot analysis of genotype × environment interaction: Proceed with caution. Crop Science, 49(5), 1564-1576.


