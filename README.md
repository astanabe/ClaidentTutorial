# ClaidentTutorial

Tutorial scripts and data for [Claident](https://github.com/astanabe/Claident).

## Prerequisites to run Claident

Claident, associated databases and programs need to be installed.
Use [ClaidentInstaller](https://github.com/astanabe/ClaidentInstaller) to install them.
In addition, [this package](https://github.com/astanabe/ClaidentTutorial/archive/main.zip) need to be downloaded and extracted to working directory.

## Prerequisites to run analyses in R

[R](https://cran.r-project.org/) and several packages need to be installed.
Install R, and then, execute the following command in R

```
library(parallel)
install.packages(c("vegan", "colorspace", "RColorBrewer", "tidyverse", "ggsci", "khroma", "picante", "bipartite", "geosphere", "foreach", "doParallel", "mpmcorrelogram", "pvclust"), repos="http://cloud.r-project.org/", dependencies=T, clean=T, Ncpus=detectCores())
```

You might need to install several programs and/or libraries such as libssl or libgdal before executing the above package installation in R.
If you are trying to run this tutorial on Ubuntu or Debian, the following commands may helps you.

```
sudo apt install libssl-dev libgdal-dev
```

## Prerequisites to learn about Claident and R

Just download [this package](https://github.com/astanabe/ClaidentTutorial/archive/main.zip) and extract to working directory.

## Learning how to run Claident

See the following shell scripts, inputs and outputs.

- [runClaident_singleend.sh](runClaident_singleend.sh)
- [runClaident_overlappedpairedend.sh](runClaident_overlappedpairedend.sh)
- [runClaident_nonoverlappedpairedend.sh](runClaident_nonoverlappedpairedend.sh)
- [runClaident_overlappedpairedend_wSTD.sh](runClaident_overlappedpairedend_wSTD.sh)

## Learning how to analyze Claident outputs in R

See the following R codes, inputs and outputs.

- [runR_overlappedpairedend.R](runR_overlappedpairedend.R)
- [runR_overlappedpairedend_wSTD.R](runR_overlappedpairedend_wSTD.R)
