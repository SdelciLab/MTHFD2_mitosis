# libraries
library(methylKit)
library(genomation)
library(reshape2)
library(ggplot2)
library(gtools)
library(tidyverse)
library(R.utils)
library(magrittr)
library(kableExtra)
library(Hmisc)
library(data.table)
library(extrafont)
library(ggbio)
library(stringr)
library(karyoploteR)
library(biovizBase)
library(lattice)
library(GenomeInfoDb)
library(BRGenomics)

# Objective: Structural analysis of Nanopore methylation data MTHFD2 WT vs KO cells

# 1. Get input data -------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

file.gff <- args[1]
file.centr <- args[2]
file.gff <- "GCF_009914755.1_T2T-CHM13v2.0_genomic.gff3.bed"
file.centr <- "GCF_009914755.1_T2T-CHM13v2.0_genomic.centr.intervals"
file.centr.bed <- "GCF_009914755.1_T2T-CHM13v2.0_genomic.centr.bed"
file.cpg_islands <- "T2T.cpg_islands.unmasked.bed"
file.genome_size <- "genome.sizes"

source("functions.R")

# 2. Pre-processing data ------------------------------------------------------------------

# column1:sample_id column2:treatment(0 or 1) column3:methylkit_file_location
design <- read.table(file = "design.txt", sep = "\t", header = TRUE)

# read centromere coords
centromeres <- read.table(text = gsub(":", "-", readLines(file.centr)),sep = "-")

# process input files, uncompress if compressed
input.list <- list()
for (i in design[,3]) {
    if ("gzfile"==summary(file(i))$class) {
        input.list <- append(input.list,
                             gunzip(filename=i, temporary=TRUE, overwrite=TRUE, remove=FALSE))
        
    } else {
        input.list <- append(input.list, i)
    }
}

# read methylation coverage for normal(5) counts, used for per base analysis
my.methraw.base <- methRead(location = input.list, sample.id=as.list(design[,1]), 
                            assembly="T2T", treatment = design[,2], context = "CpG", mincov = 5)

# read methylation coverage for low(3) counts, used for tiling and cpg islands analysis
my.methraw.low  <- methRead(location = as.list(design[,3]), sample.id=as.list(design[,1]), 
                            assembly="T2T", treatment = design[,2], context = "CpG", mincov = 3)

cpg.obj <- readFeatureFlank(file.cpg_islands, feature.flank.name=c("CpGi","shores"))
my.methraw.cpgi_and_shores <- regionCounts(my.methraw.low, GRangesList(cpg.obj$CpGi,cpg.obj$shores))
my.methraw.cpgi <- regionCounts(my.methraw.low, GRangesList(cpg.obj$CpGi))
my.methraw.shores <- regionCounts(my.methraw.low, GRangesList(cpg.obj$shores))


# 3. Quality control ------------------------------------------------------------------

my.methraw <- my.methraw.base

# plot methylation count stats per sample
pdf(file = "plots/WT.methylation.stats.pdf",  width = 8, height = 6)
getMethylationStats(my.methraw[[1]],plot=T,both.strands=F)
dev.off()
pdf(file = "plots/KO.methylation.stats.png", width = 8, height = 6)
getMethylationStats(my.methraw[[2]],plot=T,both.strands=F)
dev.off()

#plot methylation coverage per sample
pdf(file = "plots/WT.methylation.cov.png", width = 8, height = 6)
getCoverageStats(my.methraw[[1]],plot=TRUE,both.strands=FALSE,xlim=c(0.5,2.5))
dev.off()
pdf(file = "plots/KO.methylation.cov.png", width = 8, height = 6)
getCoverageStats(my.methraw[[2]],plot=TRUE,both.strands=FALSE,xlim=c(0.5,2.5))
dev.off()

# 4. Differential methylation analysis

# prepare main methylation objects
meth.obj.cpg <- methylKit::unite(my.methraw, destrand=FALSE, min.per.group=1L)
meth.obj <- meth.obj.cpg

# differential sites
myDiff.cpg <- calculateDiffMeth(meth.obj, mc.cores=1, test = "midPval")
myDiff <- myDiff.cpg

# centromeres differential methylation
myDiff.centr <- myDiff[0,]
for (i in c(seq(1,22),"X")) {
    temp <- myDiff[myDiff$chr==i & 
                       (myDiff$start>=centromeres[centromeres$V1==i,]$V2 &
                            myDiff$start<=centromeres[centromeres$V1==i,]$V3),]
    temp$chr <- paste0(i,"_centr")
    myDiff.centr <- rbind(myDiff.centr,temp)
}

# non-centromeres differential methylation
myDiff.non_centr <- myDiff[0,]
for (i in c(seq(1,22),"X")) {
    temp <- myDiff[myDiff$chr==i & 
                       (myDiff$start<centromeres[centromeres$V1==i,]$V2 |
                            myDiff$start>centromeres[centromeres$V1==i,]$V3),]
    temp$chr <- paste0(i,"")
    myDiff.non_centr <- rbind(myDiff.non_centr,temp)
}
myDiff.non_centr <- as(myDiff.non_centr,"methylDiff")

# genome-wide differential methylation
myDiff.genome <-  myDiff.non_centr
myDiff.genome$chr <- "ٴ genome"
myDiff.centr.genome <- myDiff.centr
myDiff.centr.genome$chr <- "ٴ genome_centr"

# prepare data frames
myDiff.all <- as(rbind(myDiff.non_centr, myDiff.genome, myDiff.centr, myDiff.centr.genome),"methylDiff")
myDiff.zero <-  myDiff[myDiff$meth.diff==0,]

#significant differential hyper and hypo methylated sites
myDiff15p.hyper <- getMethylDiff(myDiff,difference=15,qvalue=0.05,type="hyper")
myDiff15p.hypo <- getMethylDiff(myDiff,difference=15,qvalue=0.05,type="hypo")
myDiff15p <- getMethylDiff(myDiff,difference=15,qvalue=0.05,type="all")

#save significant methylated sites as bg
write.table(getData(myDiff15p)[,c(1,2,3,7,6)], file = "significant.methylation.all_sites.bg", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# plot significant differential hyper and hypo methylation ratio

plot_diff_meth <- function(input, title) {
    fat <- 0.50
    thin <- 0.25
    df <- input
    #remove Y
    df <- df %>% filter(chr!="Y")
    df$chr <- factor(df$chr, levels = mixedsort(levels(droplevels(df$chr))))
    df$hyper_norm <- df$percentage.of.hypermethylated / (df$percentage.of.hypermethylated+df$percentage.of.hypomethylated)
    df$hypo_norm <- df$percentage.of.hypomethylated / (df$percentage.of.hypermethylated+df$percentage.of.hypomethylated)
    df <- melt(df, id.vars = "chr",measure.vars = c(7,8))
    df$dodge <- c(rep(fat,24),rep(thin,24))
    labs <- df %>% group_by(chr) %>%summarise(value=sum(value)) %>% mutate(variable=NA)
    labs$text_lab <- paste0("(tot:",input$totals[match(labs$chr,input$chr)],
                            " ",input$number.of.hypermethylated[match(labs$chr,input$chr)],
                            " ",input$number.of.hypomethylated[match(labs$chr,input$chr)],")")
    p <- ggplot(data=df, aes(x=chr, y=value, fill=variable)) + ylim(0,1.3) +
        geom_bar(stat="identity", width = df$dodge) + coord_flip() + 
        scale_fill_discrete(name = "", labels = c("hyper", "hypo")) +
        geom_text(data=labs,aes(label=text_lab),vjust=0.2,hjust=0,size=2.5) +
        scale_fill_manual(values=c("#41BBF0","#E3781B")) +
        theme(legend.position="bottom") +
        theme_classic(base_size = 12) +
        labs(y="% (percentage)",x="Chromosome",title = title,
             subtitle = "qvalue<0.05 & methylation diff.>=25%")
    return(p)
}

#check if number hypo/hyper sites is zero and return the other for totals
checkDiffMethPerChr <- function(df) {
    for (i in 1:nrow(df)) {
        if (df[i,2]!=0) {
            df[i,6] <- 100 * df[i,2] / df[i,3]
        } else {
            if (df[i,4]!=0) {
                df[i,6] <- 100 * df[i,4] / df[i,5]
            } else {
                df[i,6] <- 0
            }
        }
    }
    return(df)
}

cairo_pdf(file = "plots/ALL.methylation.diff.all_sites.chr.pdf", width = 8, height = 6)
temp <- diffMethPerChr(myDiff.all,plot=F, keep.empty.chrom = TRUE, 
                       qvalue.cutoff=0.05, meth.cutoff=15)
temp <- temp$diffMeth.per.chr
temp$totals <- 0
temp <- checkDiffMethPerChr(temp)
plot_diff_meth(temp,"% of hyper and hypo methylated cpgs per chromosome")
dev.off()

# same plot but as pyramid

#function to plot methylation ratio as pyramid
plot_diff_meth_pyramid <- function(input, title) {
    fat <- 0.50
    thin <- 0.25
    df <- input
    #remove Y
    df <- df %>% filter(chr!="Y")
    df$chr <- factor(df$chr, levels = mixedsort(levels(droplevels(df$chr))))
    df$hyper_norm <- df$percentage.of.hypermethylated / (df$percentage.of.hypermethylated+df$percentage.of.hypomethylated)
    df$hypo_norm <- df$percentage.of.hypomethylated / (df$percentage.of.hypermethylated+df$percentage.of.hypomethylated)
    
    #add column to mark if in centromere
    df$centro <- "no"
    df$centro[grepl("centr",df$chr)] <- "yes"
    df$chr[df$centro=="yes"] <- strsplit(as.character(df$chr[df$centro=="yes"]),"_centr")
    df$chr <- droplevels(df$chr)
    df$centro <- as.factor(df$centro)
    
    #melt
    df.melt <- melt(df, id.vars = c("chr","centro"),measure.vars = c(7,8))
    
    df_pyramid <- df.melt %>% 
        mutate(value=case_when(grepl("yes",centro)~-value, TRUE~value)) %>%
        mutate(width=case_when(grepl("yes",centro)~thin, TRUE~fat))
    labs <- df_pyramid %>% group_by(chr,centro) %>%summarise(value=sum(value)) %>% mutate(variable=NA)
    labs.join.df <- left_join(labs,df,by=c("chr","centro"))
    labs$text_lab <- paste0("(tot:",labs.join.df$totals,
                            " ",labs.join.df$number.of.hypermethylated,
                            " ",labs.join.df$number.of.hypomethylated,")")
    
    p <- ggplot(df_pyramid, aes(x=value,y=chr,fill=variable)) +
        geom_bar(stat = "identity", width = df_pyramid$width) +
        scale_fill_discrete(name = "", labels = c("hyper", "hypo")) +
        geom_vline(xintercept = 0, color = "black", size=1.5) +
        scale_fill_manual(values=c("#41BBF0","#E3781B")) +
        labs(x="% (percentage)",y="Chromosome",title = title,
             subtitle = "qvalue<0.05 & methylation diff.>=25%") +
        scale_x_continuous(labels = c(1,0.5,0,0.5,1)) +
        theme_classic(base_size = 12)
    return(p)
}

cairo_pdf(file = "plots/ALL.methylation.diff.all_sites.chr.pyramid.pdf", width = 8, height = 6)
plot_diff_meth_pyramid(temp,"% of hyper and hypo methylated cpgs per chromosome")
dev.off()

# 5. Annotation analysis -----------------------------------------------------------

#read gene features
gff <- gffToGRanges(file.gff)
gene.obj <- readTranscriptFeatures(file.gff, remove.unusual = F, 
                                   up.flank = 1500, down.flank = 1000, unique.prom = F)

# parse gene features
regions.exons <- regionCounts(my.methraw,gene.obj$exons)
regions.introns <- regionCounts(my.methraw,gene.obj$introns)
regions.promoters <- regionCounts(my.methraw,gene.obj$promoters)
regions.tsss <- regionCounts(my.methraw,gene.obj$TSSes)

# annotate differential hyper and hypo methylated sites
diffAnn.hyper <- annotateWithGeneParts(as(myDiff15p.hyper,"GRanges"),gene.obj)
diffAnn.hypo <- annotateWithGeneParts(as(myDiff15p.hypo,"GRanges"),gene.obj)
diffAnn <- annotateWithGeneParts(as(myDiff15p,"GRanges"),gene.obj)

#plot annotation stats
pdf(file = "plots/ALL.methylation.diff.annot.hyper.pdf", width = 8, height = 8)
plotTargetAnnotation(diffAnn.hyper,precedence=TRUE, main="Significant differential hyper-methylation annotation")
dev.off()
pdf(file = "plots/ALL.methylation.diff.annot.hypo.pdf", width = 8, height = 8)
plotTargetAnnotation(diffAnn.hypo,precedence=TRUE, main="Significant differential hypo-methylation annotation")
dev.off()
pdf(file = "plots/ALL.methylation.diff.annot.pdf", width = 8, height = 8)
plotTargetAnnotation(diffAnn,precedence=TRUE, main="Significant differential methylation annotation")
dev.off()


# 5. Calculations and plots per cpgi -----------------------------------------------------------

# for bar plots
my.methraw.list <- list(my.methraw.cpgi_and_shores,
                        my.methraw.cpgi,
                        my.methraw.shores)
my.methraw.list.names <- list("cpgi_and_shores", "cpgi", "shores")
my.methraw.list.titles <- list("cpg islands and shores", "cpg islands", "cpg shores")

# for Karyploter
chromosomes <- read.table(file=file.genome_size,header=F,stringsAsFactors=F)
colnames(chromosomes) <- c("chromosome","length")
human.genome <- toGRanges(data.frame(chr=paste0("",chromosomes$chromosome), start=1, end=chromosomes$length))
human.genome.centromeres <- toGRanges(file.centr.bed)

for (cpgi_type in 1:3) {
    my.methraw <- my.methraw.list[[cpgi_type]]
    my.methraw.name <- my.methraw.list.names[[cpgi_type]]
    my.methraw.title <- my.methraw.list.titles[[cpgi_type]]
    
    #prepare main methylation objects
    meth.obj <- methylKit::unite(my.methraw, destrand=FALSE, min.per.group=1L)
    
    #differential sites
    myDiff <- calculateDiffMeth(meth.obj, mc.cores=1, test = "midPval")
    
    #centromeres differential methylation
    myDiff.centr <- myDiff[0,]
    for (i in c(seq(1,22),"X")) {
        temp <- myDiff[myDiff$chr==i & 
                           (myDiff$start>=centromeres[centromeres$V1==i,]$V2 &
                                myDiff$start<=centromeres[centromeres$V1==i,]$V3),]
        temp$chr <- paste0(i,"_centr")
        myDiff.centr <- rbind(myDiff.centr,temp)
    }
    #non-centromeres differential methylation
    myDiff.non_centr <- myDiff[0,]
    for (i in c(seq(1,22),"X")) {
        temp <- myDiff[myDiff$chr==i & 
                           (myDiff$start<centromeres[centromeres$V1==i,]$V2 |
                                myDiff$start>centromeres[centromeres$V1==i,]$V3),]
        temp$chr <- paste0(i,"")
        myDiff.non_centr <- rbind(myDiff.non_centr,temp)
    }
    myDiff.non_centr <- as(myDiff.non_centr,"methylDiff")
    
    #genome-wide differential methylation
    myDiff.genome <-  myDiff.non_centr
    myDiff.genome$chr <- "ٴ genome"
    myDiff.centr.genome <- myDiff.centr
    myDiff.centr.genome$chr <- "ٴ genome_centr"
    
    #prepare data frames
    myDiff.all <- as(rbind(myDiff.non_centr, myDiff.genome, myDiff.centr, myDiff.centr.genome),"methylDiff")
    myDiff.zero <-  myDiff[myDiff$meth.diff==0,]
    myDiff25p <- getMethylDiff(myDiff,difference=15,qvalue=0.05,type="all")
    
    #save significant methylated sites as bg
    write.table(getData(myDiff25p)[,c(1,2,3,7,6)], 
                file =paste0("significant.methylation.",my.methraw.name,".bg"), 
                quote = F, sep = "\t", row.names = F, col.names = F)
    
    #plot significant differential hyper and hypo methylation cpgi
    cairo_pdf(file = paste0("plots/ALL.methylation.diff.",my.methraw.name,".chr.pdf"), 
              width = 8, height = 6)
    temp <- diffMethPerChr(myDiff.all,plot=F, keep.empty.chrom = TRUE,
                           qvalue.cutoff=0.05, meth.cutoff=15)
    temp <- temp$diffMeth.per.chr
    temp$totals <- 0
    temp <- checkDiffMethPerChr(temp)
    print({
        p <- plot_diff_meth(temp,paste0("% of hyper and hypo methylated ",my.methraw.title," per chromosome"))
    })
    dev.off()
    
    #same plot but as pyramid
    cairo_pdf(file = paste0("plots/ALL.methylation.diff.",my.methraw.name,".chr.pyramid.pdf"), 
              width = 8, height = 6)
    print({
        p <- plot_diff_meth_pyramid(temp,paste0("% of hyper and hypo methylated ",my.methraw.title," per chromosome"))
    }) 
    dev.off()
    
    #karyploter
    meth.cpgi_ranges <- toGRanges(paste0("significant.methylation.",my.methraw.name,".bg"))
    meth.cpgi_ranges.hyper <- meth.cpgi_ranges[meth.cpgi_ranges$V4>0]
    meth.cpgi_ranges.hypo  <- meth.cpgi_ranges[meth.cpgi_ranges$V4<0]
    
    cairo_pdf(file = paste0("plots/significant.methylation.",my.methraw.name,".karyoplot.pdf"), 
              width = 8, height = 6)
    kp <- plotKaryotype(genome = human.genome,cex=0.7,plot.type=1,
                        main=my.methraw.title)
    kpRect(kp,data=human.genome.centromeres,y0=0, y1=0.25, r0=-0.4, col="#AAAAFF")
    kpBars(kp,data=meth.cpgi_ranges.hyper, y1=meth.cpgi_ranges.hyper$V4/100*0.5, r0=-0.05, col="#41BBF0", border="#41BBF0")
    kpBars(kp,data=meth.cpgi_ranges.hypo, y1=meth.cpgi_ranges.hypo$V4/100*0.5, r0 = -0.4, col="#E3781B", border="#E3781B")
    dev.off()
}
