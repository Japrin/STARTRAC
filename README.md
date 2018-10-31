[![Travis CI Build Status](https://travis-ci.org/Japrin/STARTRAC.svg?branch=master)](https://travis-ci.org/Japrin/STARTRAC)

# STARTRAC
STARTRAC(Single T-cell Analysis by Rna-seq and Tcr TRACking).

# Installation

To install this package, simply:
```
install.packages("devtools")
devtools::install_github("Japrin/STARTRAC")
```

# Example
Read input data:
```
dat.file <- system.file("extdata/example.cloneDat.Zhang2018.txt",package = "Startrac")
in.dat <- read.table(dat.file,stringsAsFactors = F,head=T)
```

run the STARTRAC pipeline:
```
out <- Startrac.run(in.dat, proj="CRC", cores=NULL,verbose=F)
```

More information can be found in the [vignettes](http://htmlpreview.github.io/?https://github.com/Japrin/STARTRAC/blob/master/vignettes/startrac.html).
If you use this package, please cite [this Nature paper](https://www.nature.com/articles/s41586-018-0694-x). Thank you! 

