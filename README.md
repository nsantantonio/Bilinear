# Bilinear

Fit bilinear models using AMMI (Additive Main effects Multiplicative Interaction) or GGE/SREG (Genotype and Genotype by Environment /Sites REGression) and make 2D and 3D biplots.  Also tests for the number of significant dimensions of GxE using a parameteric bootstrap test method suggested by Forkman and Piepho (2014). An additional test for significant multiplicative terms, the F<sub>R</sub> test from Piepho (1995), has been also been implemented in version 1.8.5.

See 'Bilinear examples.R' for a tutorial of use.  ```Bilinear vx.x.R``` contains the relevant functions, namely ```Bilinear()``` to fit model and ```BBplot()``` to plot 2D and 3D biplots. With version 1.8.5, an additional funciton ```AMMIplot()```, can be used to produce linear and winner plots such as those in Gauch & Zobel (1997).

### Plant Breeders Disclaimer
It is this author's opinion that GGE type models are more difficult to interpret than AMMI type models, as the genotypic main effect will be present in one or more of the remaining dimensions of the residual matrix, but not necessarily the first (largest) dimension (e.g. when Var(GxE) > Var(G)).  Therefore, I recommend using AMMI type models, as both genotype and environment main effects are removed from the residuals before GxE effects are assessed.  However, GGE type models are provided here as a freely available resource for those who want to use them. The ```decorateGGE``` option of BBplot will draw the mega-environment delineations suggested by Yan et al. (2000).


## Known issues to be addressed

* Need balanced data.  At this time, the program requires balanced data across genotypes and environments (i.e. all genotypes observed in all environments).  If just one or two cells are missing, you could impute the genotype effect + environment effect (i.e. no GxE) for that cell. The program will run with unequal replication within each location, but each genotype must be observed at least once in each environment, and unequal replication could result in erroneous estimates (the program should print a warning if there is unequal replication).  Eventually an EM algorithm might be implemented to account for unbalanced data as suggested by Gauch and Zorbel (1990). This is planned to be implemented soon.

### Eventually R package to CRAN ?:
If this author finds time to produce a more polished, flexible and tested program, it may be submitted to CRAN as an R package. Until then, please follow me here on github and send me a message, particularly if you decide to publish any results produced from this (NOT because I want authorship, I just want to know if and how it is being used, particularly for selections)

### Soy data
The ```soy.Rdata``` included here is from Zobel, Wright and Gauch (1988) and is subset from the ```gauch.soy``` dataset included in the [```agridat``` package](https://github.com/kwstat/agridat). This was included to show the equivalence between using within environment means and raw data. 

## References 
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
- Yang, R. C., Crossa, J., Cornelius, P. L., & Burgueño, J. (2009). Biplot analysis of genotype× environment interaction: Proceed with caution. Crop Science, 49(5), 1564-1576.





