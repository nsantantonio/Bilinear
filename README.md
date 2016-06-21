# Bilinear
Fit bilinear models using AMMI or GGE and make 2D and 3D biplots

See 'Bilinear examples.R' for a tutorial of use.  'Bilinear vx.x.R' contains the relevant functions.

Known issues to be addressed:

- Need balanced data.  At this time, the program requires balanced data across genotypes and environments (i.e. all genotypes observed in all environments).  If just one or two cells are missing, you could impute the genotype effect + environment effect (i.e. no GxE). The program will run with unequal replication within each location, but each genotype must be observed at least once in each environment, and unequal replication could result in erroneous estimates (the program should print a warning if there is unequal replication).  

- SVD sign ambiguity. The sign of Genotype and Environment scores is not necessarily the sign of the effect.  Still working on a resolution.
 
- Flexability to change point types, colors, etc.
