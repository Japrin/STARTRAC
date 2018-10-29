# STARTRAC
STARTRAC(Single T-cell Analysis by Rna-seq and Tcr TRACking)

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

More information can be found in the [vignettes](http://htmlpreview.github.io/?https://github.com/Japrin/STARTRAC/blob/master/vignettes/startrac.html)

