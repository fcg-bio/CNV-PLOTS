#!/usr/bin/Rscript
#/opt/well/R/R-3.0.2/bin/Rscript

cran.packages <- c(
	"RColorBrewer"
	,"ggplot2"
	)
biocLite.packages <- c("Biostrings", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges", "ggbio", "biovizBase", "gridExtra")

pkgTest <- function(x, y)
{
	if (!require(x,character.only = TRUE)) 
	{
	cat(paste("Installing package ", x ,"\n\n", sep=""))	
	if(y=="CRAN") install.packages(x,dep=TRUE,repos="http://www.stats.bris.ac.uk/R")
	if(y=="biocLite") source("http://bioconductor.org/biocLite.R")
	if(y=="biocLite") biocLite(x, ask=F)
        if(!require(x,character.only = TRUE)) stop(x, "Package not found")
	}
}

lapply(cran.packages, pkgTest, y="CRAN")
lapply(biocLite.packages, pkgTest, y="biocLite")

all.packages <- c(cran.packages,biocLite.packages)
check.pckg <- all.packages %in% installed.packages()
if(! any(check.pckg==FALSE)) print("Congratulations: All the packages are installed!!!")
if(any(check.pckg==FALSE)) print(paste("Packages not installed:", check.pckg[!check.pckg], sep=""))

