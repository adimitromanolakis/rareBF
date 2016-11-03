# JLS

R package for  the Joint Location Scale (JLS) test (Soave et al. 2015) to simultaneously test for mean and variance differences between groups. The JLS test uses Fisher's combined p-value method to combine evidence from the individual locaiton (regression t-test) and scale (Levene's test of homogeneity of variances) tests.  Author: David Soave (david.soave@mail.utoronto.ca).

See the reference manual (JLS.pdf) for package details.

##Installation Instructions

Users can install the most recent version of the JLS package using install_github() in the devtools package.

```
install.packages("devtools")
devtools::install_github("dsoave/JLS")
```
```
library(JLS)
?JLS_test
```





