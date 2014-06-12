#!/usr/bin/Rscript
#/opt/well/R/R-3.0.2/bin/Rscript


# Usage
{
  usage <- "
  cnv-plot.R : R script to create plots of cnv distribution across genome in multiple samples. Planned plots: 
    1) Barplot showing gain and loses, with bars representing the number of samples having a CNV change.
    2) Tile plot detailing the cnv changes by sample.
    3) Join plot (1+2) using the new multiplot function
  Usage: Rscript cnv-plot.R [options] <input file>
  Example : cnv-plot.R prefix=cnv-plot file.tsv

  Input format: Tab separated file with the column names (the order is not important)
    1 Chromosome REQ
    2 Start REQ
    3 End REQ
    4 Sample REQ : Sample names
    5 Type REQ : CNV type. Possible values (cap or small caps) GAIN, DELETION, LOSS 
    6 CopyNumber OPT : Copy number in the region
    7 LOH OPT : Loss of heterozygosity in the region
		
  Options:
    1 bin : Bin size for windows. Default = 1000000
    2 main : Main for plots. Default = Empty
    3 bsgenome : BSgenome library. Default = BSgenome.Hsapiens.UCSC.hg19
    4 prefix : Prefix for output files. Default = cnv-plot
	"
}


# Reading and Processing argments
{
# Reading arguments
argsRaw<-commandArgs(trailingOnly = T)
if (length(argsRaw) == 0){ cat(usage,"\n"); stop("File not provided. Check Usage")}

# Creating argument list
arguments <- list (file=NULL, main="", bin=1000000, BSgenome="BSgenome.Hsapiens.UCSC.hg19", prefix = "cnv-plot")

# Get input file name
arguments$file = argsRaw[length(argsRaw)]

if(!file.exists(arguments$file)) {at(usage,"\n");stop("Input file doesn't exists")}

# Read the other arguments
for(i in 1:(length(argsRaw)-1)){
  argsRaw2<-unlist(strsplit(argsRaw[i],"="))
  arguments[argsRaw2[1]]<-argsRaw2[2]
}
}

# Reading and processing input files
{
  dat <- read.delim(arguments$file, stringsAsFactors = F)

  # Check Chromosome Column Format
  dat$Chromosome <- gsub("chr","",dat$Chromosome)
  
  # Check Type Column Format
  dat$Type <- tolower(dat$Type)
  errorFormat <- which(!dat$Type %in% c("gain","deletion","loss"))
  if(length(errorFormat)>0) {cat(usage,"\n"); stop("Please provide a correct values for type column")}
  dat$Type[dat$Type == "loss"] <- "deletion"
  
  # Check Other Column Formats
  if(!class(dat[,"Start"])=="integer") {cat(usage,"\n"); stop("Please provide a correct Start input format (integer)")}
  if(!class(dat[,"End"])=="integer") {cat(usage,"\n"); stop("Please provide a correct End input format (integer)")}
  
  
}


CNV.sample.tile.plot <- function(data, bin=1000000, BSgenome="BSgenome.Hsapiens.UCSC.hg19", main="") {
  # Description
  # From a list of CNV in multiple samples it performs a tile plot with gains and loses
  
  # Usage
  # CNV.sample.count.plot(data, bin, main="", BSgenome="BSgenome.Hsapiens.UCSC.hg19")
  
  # Arguments
  # data Format: Required columns Chromosome, Start, End, Sample, Type

  # Check default parameters
  if(!exists("bin")) bin <- 1000000
  if(!exists("BSgenome")) BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
  if(!exists("main")) main <- ""
  
  # Loading Libraries
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(BSgenome, character.only = T))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggbio))
  suppressPackageStartupMessages(library(gridExtra))
  
  # Arguments format
  bin <- as.numeric(bin)
  
  # Get Chromosomes lengths
  chrlengths <- seqlengths(get(BSgenome))
  names(chrlengths) <- gsub("chr", "", names(chrlengths))
  chrNumeric <- grep("^[0-9]+$", names(chrlengths), value = T) # Get only numeric chr
  chrNonNumeric <- grep("^[0-9]+$", names(chrlengths), value = T, invert = T)
  chrNonNumeric <- chrNonNumeric[chrNonNumeric %in%  data$Chromosome]
  chrLevels <- c(chrNumeric, chrNonNumeric)
  chrlengths <- chrlengths[chrLevels]
  
  # Convert coordinates by bin size
  data$Start <- ceiling(as.integer(data$Start)/bin)
  data$End <- ceiling(as.integer(data$End)/bin)
  chrlengths <- ceiling(chrlengths/bin)
  
  # Creating positions matrix
  if(exists("posM")){ rm(posM)}
  for (chr in names(chrlengths)) {
    for (i in 1:chrlengths[chr]){
      if(exists("posM")) {
        posM <- rbind(posM, c(chr,i))
      } else {
        posM <- c(chr,i)
      }
    }  
  }
  row.names(posM) <- NULL
  posM <- data.frame(posM)
  names(posM) <- c("chr","bin")
  posM$id <- paste(posM$chr,posM$bin, sep="_")
  posM$pos <- 1:dim(posM)[1]
  posM <- posM[,c("id","pos")]

  # Add position to Data
  data$StartBin <- paste(data$Chromosome, data$Start, sep="_")
  data$EndBin <- paste(data$Chromosome, data$End, sep="_")
  data <- merge(data, posM, by.x="StartBin", by.y="id")
  data$StartPos <- data$pos ; data$pos <- NULL
  data <- merge(data, posM, by.x="EndBin", by.y="id")
  data$EndPos <- data$pos ; data$pos <- NULL
  
  # Remove unnecessary columns that can produce duplicates
  data <- unique(data[,c("Chromosome","Start","End","Sample","Type","StartPos","EndPos")])
  
  # Dealing with deletions spanning more than 1 bin
  y1 <- data[data$StartPos == data$EndPos,]
  y2 <- data[data$StartPos != data$EndPos,]
  
  ly <- apply(y2,1,function(yl){
    npos <- as.numeric(yl["StartPos"]):as.numeric(yl["EndPos"])
    if(exists("res")) rm(res)
    for (i in npos){
      out <- c(yl["Chromosome"],yl["Start"],yl["End"],yl["Sample"],yl["Type"],i,i)
      if(exists("res")) res <- rbind(res, out)
      if(!exists("res")) res <- out
    }
    res
  })
  
  lymat <- do.call(rbind, ly)
  dimnames(lymat)[[2]] <- names(y1)
  dimnames(lymat)[[1]] <- 1:dim(lymat)[1]
  lymat <- data.frame(lymat)
  row.names(lymat) <- NULL
  row.names(y1) <- NULL
  data <- rbind(y1,lymat)
  data$StartPos <- as.numeric(data$StartPos)
  
  # Plot
  sl <- chrlengths
  csl <- cumsum(chrlengths)
  breaks <- csl - (sl/2)
  data$TypeBinary[data$Type=="gain"] <- 1
  data$TypeBinary[data$Type=="deletion"] <- -1
  xlim <- c(1, max(posM$pos))

  p <- ggplot(data, aes(x=StartPos, y=Sample, fill=factor(Type)))+
    geom_tile()+
    scale_x_continuous(name = "Chromosome", breaks=breaks, minor_breaks = csl[-length(csl)]+0.5, limits=xlim, expand=c(0,0))+
    scale_fill_manual(values=c("#C93312","#899DA4"))+
    theme(panel.background=element_blank())+
    theme(legend.position="none") +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())+
    theme(panel.grid.minor.x=element_line(colour='red',linetype='dashed')) + 
    ggtitle(main)
  
  p

}


CNV.sample.count <- function(data, bin=1000000, BSgenome="BSgenome.Hsapiens.UCSC.hg19", main="") {
  # Description
  # From a list of CNV in multiple samples it calculates the number of gains and deletions in each bin and plot the summary in Barplot. If a sample has more than one CNV coordinate in the bin, it counts only once
  
  # Usage
  # CNV.sample.count.plot(data, bin, main="", BSgenome="BSgenome.Hsapiens.UCSC.hg19")
  
  # Arguments
  # data Format: Required columns Chromosome, Start, End, Sample, Type
  
  # Loading Libraries
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(BSgenome, character.only = T))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggbio))
  suppressPackageStartupMessages(library(biovizBase))
  suppressPackageStartupMessages(library(gridExtra))
  
  # Arguments format
  bin <- as.numeric(bin)
  
  # Get only necessary columns from the data
  data <- data[,c("Chromosome","Start","End","Sample","Type")]
  
  # Get Chromosomes lengths
  chrlengths <- seqlengths(get(BSgenome))
  names(chrlengths) <- gsub("chr", "", names(chrlengths))
  chrNumeric <- grep("^[0-9]+$", names(chrlengths), value = T) # Get only numeric chr
  chrNonNumeric <- grep("^[0-9]+$", names(chrlengths), value = T, invert = T)
  chrNonNumeric <- chrNonNumeric[chrNonNumeric %in%  data$Chromosome]
  chrLevels <- c(chrNumeric, chrNonNumeric)
  chrlengths <- chrlengths[chrLevels]
  
  # Correction for missing Chromosome data 
  # fake entries that must be set to coverage = 0 in res data.frame
  # All combinations Type - Chromosome
  cnvDelete <- list(c(),c())
  names(cnvDelete) <- levels(factor(data$Type))
  for (typeCnv in levels(factor(data$Type)) ) {
    sdata <- subset(data, Type==typeCnv)
    emptyChr <- chrLevels[!chrLevels %in% sdata$Chromosome]
    for (eChr in emptyChr) {
      cnvDelete[[typeCnv]] <- c(cnvDelete[[typeCnv]],eChr)
      data[dim(data)[1]+1,c("Chromosome","Start","End","Type")] <- c(eChr,1,1,typeCnv)
    }
  }
  
  # Convert coordinates by bin size
  data$Start <- ceiling(as.integer(data$Start)/bin)
  data$End <- ceiling(as.integer(data$End)/bin)
  chrlengths <- ceiling(chrlengths/bin)
  
  # Get only unique values from data. If a sample has more than one CNV coordinate in the bin, it counts only once
  data <- unique(data)
  
  # Conver to Factor
  data$Chromosome <- factor(data$Chromosome, chrLevels)
  data$Type <-  factor(data$Type, levels=c("gain","deletion"))
  
  # Create GenomicRanges Object
  gr <- GRanges(seqnames = data$Chromosome, IRanges(start = data$Start, end = data$End), Type = data$Type)
  seqlengths(gr) <- chrlengths
  facets <- as.formula(Type~seqnames)
  grl <- splitByFacets(gr, facets)
  
  # Sort grl by chromosome and type
  z <- outer(levels(data$Type), levels(data$Chromosome), paste, sep=".")
  z <- c(z[1,],z[2,])
  z <- z[z %in% paste(data$Type, data$Chromosome, sep=".")]
  grl <- grl[z]
  
  # Create Coverage Data Frame
  xlim <- c(min(start(ranges(gr))),max(end(ranges(gr))))
  lst <- lapply(grl, function(dt){
    vals <- coverage(keepSeqlevels(dt, unique(as.character(seqnames(dt)))))
    seqName <- as.character(seqnames(dt)@values)
    if(any(is.na(seqlengths(dt)))){
      seqs <- xlim[1]:max(end(dt))
      vals <- vals[[1]][seqs]
      vals <- as.numeric(vals)                           
      vals <- c(vals, rep(0, xlim[2]-max(end(dt))))
      seqs <- xlim[1]:xlim[2]
    }else{
      seqs <- 1:seqlengths(dt)[seqName]
      vals <- vals[[1]][seqs]
      vals <- as.numeric(vals)                           
    }
    if(length(unique(values(dt)$.id.name))){                
      res <- data.frame(coverage = vals, seqs = seqs,
                        seqnames =
                          as.character(seqnames(dt))[1],
                        .id.name = unique(values(dt)$.id.name))
    }else{
      res <- data.frame(coverage = vals, seqs = seqs,
                        seqnames = as.character(seqnames(dt))[1])      
    }
    res
  })
  res <- do.call(rbind, lst)
  lres <- lapply(row.names(res), function(x){unlist(strsplit(x, "\\."))[1:2]})
  lrmat <- do.call(rbind, lres)
  dimnames(lrmat) <- list(c(),c("Type","Chromosome"))
  res <- cbind(res,lrmat)
  res$Type <- relevel(res$Type, ref="gain")
  res$Position <- NA
  res$Position[res$Type=="gain"] <- 1:dim(res[res$Type=="gain",])[1]
  res$Position[res$Type=="deletion"] <- 1:dim(res[res$Type=="gain",])[1]
  res$Chromosome <- factor(res$Chromosome, chrLevels)
  res$coverage[res$Type=="gaib" & res$Chromosome %in% cnvDelete$gain] <-0
  res$coverage[res$Type=="deletion" & res$Chromosome %in% cnvDelete$deletion] <-0
  res$coverage[res$Type=="deletion"] <- -res$coverage[res$Type=="deletion"]
  res$seq_type <-paste(res$Type, res$Chromosome, sep="_")
  res$seq_type <- factor(res$seq_type, levels=unique(res$seq_type))
  
  # Plot
  sl <- seqlengths(grl)
  csl <- cumsum(seqlengths(grl))
  breaks <- csl - (sl/2)
  cols <- c(rep(c("#899DA4","#78B7C5"),nlevels(res$Chromosome)/2), rep(c("#C93312","#F21A00"),nlevels(res$Chromosome)/2))
  names(cols) <- levels(res$seq_type)
  ylimits <- c(-max(abs(res$coverage)), max(abs(res$coverage)))
  
  p <- ggplot(data = res, aes(x=Position, width=1)) + 
    scale_x_continuous(name = "Chromosome", breaks = breaks, minor_breaks = csl[-length(csl)]+0.5,  expand=c(0,0)) + 
    scale_y_continuous(limits=ylimits) +
    geom_bar(data=subset(res, Type=="gain"), aes(y=coverage, fill=seq_type), stat = "identity", position = "identity") + 
    geom_bar(data=subset(res, Type=="deletion"), aes(y=coverage, fill=seq_type), stat = "identity", position = "identity") + 
    geom_hline(aes(yintercept=0), colour="#FFFFFF", size = 1) + 
    theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks = element_blank()) +
    theme(panel.grid.minor.x=element_line(colour='red',linetype='dashed'))+
    theme(panel.background=element_blank()) +
    theme(legend.position="none") + 
    ylab(paste("Deletion","Gain", sep="      ")) +
    scale_fill_manual(values = cols) + 
    ggtitle(main)
  # p + annotate("text", x = breaks, y = 0, label = names(breaks))  
  
  # Return plot and table
  return(list("plot" = p, "table" = res))
  
}


p1 <- CNV.sample.count(dat, bin=arguments$bin, main=arguments$main, BSgenome=arguments$BSgenome)


p2 <- CNV.sample.tile.plot(data=dat, bin=arguments$bin, main=arguments$main, BSgenome=arguments$BSgenome)

write.table(p1$table, file = paste(arguments$prefix,"sum.txt", sep="."), row.names=F, quote=F, sep="\t")

pdfName <- paste(arguments$prefix,"sum.pdf", sep=".")
pdf(pdfName, width= 19, height= 8, title=arguments$main)
print(p1$plot)
dev.off()

pdfName <- paste(arguments$prefix,"bySample.pdf", sep=".")
pdf(pdfName, width= 19, height= 8, title=arguments$main)
print(p2)
dev.off()

pdfName <- paste(arguments$prefix,"combined.pdf", sep=".")
pdf(pdfName, width= 19, height= 16, title=arguments$main)
grid.arrange(p1$plot, p2, ncol=1, nrow=2, heights=c(1, 4))
dev.off()
