# ClaidentTutorial

Tutorial scripts and data for Claident.

## Prerequisites to run Claident

Claident, associated databases and programs need to be installed.
Use [ClaidentInstaller](https://github.com/astanabe/ClaidentInstaller) to install them.
In addition, [this package](https://github.com/astanabe/ClaidentInstaller/archive/master.zip) need to be downloaded and extracted to working directory.

## Prerequisites to run analysis in R

[R](https://cran.r-project.org/) and several packages need to be installed.
Install R, and then, execute the following command in R

```
library(parallel)
install.packages(c("vegan", "colorspace", "RColorBrewer", "tidyverse", "ggsci", "khroma", "picante", "bipartite", "geosphere", "foreach", "doParallel", "mpmcorrelogram"), repos = "http://cloud.r-project.org/", dependencies=T, clean=T, Ncpus=detectCores())
```

## Prerequisites to learn about Claident and R

Just download [this package](https://github.com/astanabe/ClaidentInstaller/archive/master.zip) and extract to working directory.
