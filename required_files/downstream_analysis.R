### Secondary analysis of ATACseq samples
### Change working directory into the directory with the "final.bam" files
library(tidyverse)
library(RCurl)
library(ggsci)
library(Rsamtools)
library(GenomicAlignments)
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/plot_cleanup.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/convert_to_human_readable_numbers.R", ssl.verifypeer = FALSE)))


### Created list of .bam files and sample info ##################################################################
bamFiles <- dir(path = "..", pattern = "final.bam$", full.names = TRUE)
sampleList <- read_tsv("sample_list.txt")  ### create a tab separated list of bam file names, sample names and group names

### Function to plot read numbers mapped to each chromosome ######################################################
plotChromReads <- function(bamList = bamFiles){

  chrom <- tibble(Chromosome = idxstatsBam(bamList[1])[[1]], Length = idxstatsBam(bamList[1])[[2]])
  for (i in 1:dim(sampleList)[1]){
    chrom <- cbind(chrom, idxstatsBam(bamList[i])[3])
    colnames(chrom)[i+2] <- sampleList[i,2]
  }
  chrom$Chromosome <- gsub("_C_albicans_SC5314", "", chrom$Chromosome)
  chrom$Chromosome <- gsub("Ca22", "", chrom$Chromosome)
  ChromReadPlot <- chrom %>% filter(Chromosome != "chrM") %>% gather(key = "Sample", value = "Reads", -Chromosome, -Length) %>%
    ggplot(aes(x = Chromosome, y = Reads, fill = Sample, group = Sample)) + 
    geom_bar(stat = "identity", position = "dodge") + coord_flip() +
    ylab("No. of mapped reads") +
    cleanup + scale_fill_npg()
  return(ChromReadPlot)
}

plotChromReads()
pdf(file = "chromReadsPlot.pdf")
plotChromReads()
dev.off()

png(file = "chromReadsPlot.png")
plotChromReads()
dev.off()


### Function to extract the read length distributions of all samples ##############################################

extractFragmentSize <- function(bamList = bamFiles, sampList = sampleList){
  
  insertSizesAll <- read.table(text = "", colClasses = c("numeric", "numeric", "factor"), 
                             col.names = c("insertSizes", "Freq", "sample"))
  for (i in 1:dim(sampList)[1]){
    atacReads <- readGAlignmentPairs(bamList[i], param = ScanBamParam(what = c("qname", "isize")))
    # param = ScanBamParam(mapqFilter = 1, flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE) can be used as additional paramters for filtering
    atacReads_read1 <- GenomicAlignments::first(atacReads)
    rm(atacReads)
    gc(verbose = FALSE)
    sname <- sampList[[i,2]]
    insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
    rm(atacReads_read1)
    gc(verbose = FALSE)
    # procude a dataframe with isizes and add a column with the sample name, append every round in the loop to this dataframe
    insertSizesDF <- data.frame(table(insertSizes))
    insertSizesDF$sample <- sname
    insertSizesAll <- rbind(insertSizesAll, insertSizesDF)
  }
  return(insertSizesAll)
}

insertSizesAll <- extractFragmentSize()

### Function to plot the read length distributions of all samples #####################################################

plotFragmentSize <- function(sampList = sampleList[[2]], isizes = insertSizesAll){  
  ### character vector of sample names from sample list can be used as sampList argument to chose specific samples
  
  fragLenPlot <- isizes %>% 
    mutate(InsertSize = as.numeric(as.vector(insertSizes)), Count = as.numeric(as.vector(Freq)), 
         Sample = as.factor(as.vector(sample))) %>% 
    .[.$Sample %in% sampList,] %>%
    ggplot(aes(x = InsertSize, y = Count, color = Sample)) + 
    geom_line() + cleanup + scale_color_npg() + theme(legend.key = element_rect(fill = NA))
  return(fragLenPlot)
}

fragLenPlot <- plotFragmentSize()
fragLenPlot

pdf(file = "fragLengthPlot.pdf", height = 5, width = 7.5)
plotFragmentSize()
dev.off()

png(file = "fragLengthPlot.png", height = 1440, width = 2160, res = 300)
plotFragmentSize()
dev.off()


### Function to calculate number of reads within peaks for each sample ################################################
library(rtracklayer)
library(ChIPQC)

Ca_feat <- import("../required_files/C_albicans_SC5314_A22_current_features_haploid.gff")
CaAnno <- list(version="custom", Ca_features=Ca_feat)
openRegionPeaks <- dir(path = "../macs2_peaks_calling/", pattern = "narrowPeak$", full.names = TRUE)
nucFreeFiles <- dir(path = "..", pattern = "bam.nucfree.bam$", full.names = TRUE)
sampleList_wo_contr <- sampleList[sampleList$group != "gDNA", ]

calcReadsInPeaks <- function(bamFiles = nucFreeFiles, Peaks = openRegionPeaks){
  
  peakCallResult <- read.table(text = "", colClasses = c("numeric", "numeric", "character"), 
                             col.names = c("Reads", "RiP.", "Peaks", "Sample"))

  for (i in 1:length(bamFiles)) {
    qcRes <- ChIPQCsample(bamFiles[i], peaks = Peaks[i], annotation = CaAnno)
    res <- QCmetrics(qcRes) %>% t %>% data.frame %>% 
        dplyr:::select(Reads, starts_with(c("RiP")))
    res$Peaks <- length(qcRes$Counts)
    res$Sample <- sampleList_wo_contr$sample[i]
    peakCallResult <- rbind(peakCallResult, res)
  }
  return(peakCallResult)
}

peakCallResult <- calcReadsInPeaks()

RiPplot <- peakCallResult %>% mutate(RiP. = RiP.*100) %>% 
  gather(key = "Data", value = "Peaks", RiP.:Peaks) %>%
  ggplot() + geom_col(aes(x = Sample, y = Peaks, fill = Data), position = "dodge") +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "% of reads in peaks"), name = "No. of peaks") +
  labs(x = "") + coord_flip() + cleanup + scale_fill_npg(name = "", labels = c("No. of peaks", "% of reads in peaks"))

RiPplot

pdf(file = "ReadsInPeaksPlot.pdf", height = 4)
RiPplot
dev.off()

png(file = "ReadsInPeaksPlot.png", res = 300, width = 1920, height = 960)
RiPplot
dev.off()


### Function to annotate peaks called by MACS2 ######################################################################
library(rtracklayer)
library(ChIPseeker)
library(GenomicFeatures)

mergedPeaks <- import("../macs2_peaks_calling/merged_all_peaks_annot.bed")
gff_file <- "../required_files/C_albicans_SC5314_A22_current_features_haploid.gff"

annotatePeaks <- function(peaks = mergedPeaks){
  TxDb.Calbicans.A22.CGD.features <- makeTxDbFromGFF(gff_file, dataSource = "CGD", organism = "Candida albicans",
                                                     circ_seqs = "Ca22chrM_C_albicans_SC5314")
  anno <- annotatePeak(mergedPeaks, TxDb = TxDb.Calbicans.A22.CGD.features)
  return(anno)
}

mergedPeaksAnno <- annotatePeaks()

plotAnnoPie(mergedPeaksAnno)
plotAnnoBar(mergedPeaksAnno)
upsetplot(mergedPeaksAnno)

pdf(file = "Peak_annotation_plots.pdf", height = 4)
upsetplot(mergedPeaksAnno)
plotAnnoBar(mergedPeaksAnno)
plotAnnoPie(mergedPeaksAnno)
dev.off()

png(file = "upsetplot_peak_annotation.png", res = 300, width = 1920, height = 1440)
upsetplot(mergedPeaksAnno)
dev.off()
png(file = "barplot_peak_annotation.png", res = 300, width = 1920, height = 720)
plotAnnoBar(mergedPeaksAnno)
dev.off()
png(file = "pieplot_peak_annotation.png", res = 300, width = 1920, height = 1440)
plotAnnoPie(mergedPeaksAnno)
dev.off()


### Venn diagram creation to visualize the overlaps in peak calling between samples #############

convertPeakBed <- function(mPeaks = mergedPeaks, samplList = sampleList_wo_contr){

  for (i in 1:length(samplList[[2]])){
    mPeaks$x <- ifelse(grepl(gsub(samplList[i,1], 
                                          pattern = ".bam.fq.final.bam", 
                                          replacement = ""), mPeaks$name), 1, 0)
    sname <- samplList[[i,2]]
    colnames(mcols(mPeaks))[i + 1] <- sname
  }
  mPeaks$name <- NULL
  return(mPeaks)
}

mergedPeaksSample <- convertPeakBed()

library(limma)

pdf(file = "VennDiagrams_Peaks.pdf", height = 5)
as.data.frame(elementMetadata(mergedPeaksSample)) %>% dplyr::select(starts_with("YPD")) %>% 
    vennDiagram(main = "Overlap for YPD nucleosome free regions", lwd = 1.5, cex=1.2)
as.data.frame(elementMetadata(mergedPeaksSample)) %>% dplyr::select(starts_with("H2O2")) %>% 
    vennDiagram(main = "Overlap for H2O2 nucleosome free regions", lwd = 1.5, cex = 1.2)
dev.off()

#as.data.frame(elementMetadata(mergedPeaksSample)) %>% dplyr::select(ends_with("_2")) %>%  ### use this for within replicate comparison
 #   vennDiagram(main = "Overlap for replicate 2 nucleosome free regions", lwd = 1.5, cex = 1.2)


PCApeaks <- as.data.frame(elementMetadata(mergedPeaksSample)) %>%
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("_\\d", "", Samples))
  
adj <- ifelse(PCApeaks$PC1 > mean(PCApeaks$PC1), 1.2, -0.2)

PCAplot_peaks <- ggplot(data = PCApeaks, aes(x = PC1, y = PC2, colour = Group)) + 
  theme(legend.key = element_rect(fill = NA)) +
  geom_text(aes(label = Samples), hjust = adj) +
  geom_point(size = 5) + cleanup + scale_color_npg()

PCAplot_peaks

pdf(file = "PCAplot_peaks.pdf", height = 5)
PCAplot_peaks
dev.off()

png(filename = "PCAplot_peaks.png", res = 300, height = 1440, width = 1920)
PCAplot_peaks
dev.off()





#peakFiles <- dir(path = "../macs2_peaks_calling/", pattern = "narrowPeak$", full.names = TRUE)
#nucFreeFiles <- dir(path = "..", pattern = "bam.nucfree.bam$", full.names = TRUE)
#peakList <- sampleList[sampleList$group != "gDNA", ]
#peakList <- tibble(SampleID = peakList[[2]], Group = peakList[[3]], bamReads = nucFreeFiles, Peaks = peakFiles)
#nucFreeExp <- ChIPQC(as.data.frame(peakList), annotation = CaAnno) ### needs too much memory, not possible

#library(DT)
#library(GenomicFeatures)
#Ca_genes <- Ca_feat[Ca_feat$type == "gene",]
#CaATGs <- resize(Ca_genes, fix = "start", 1)
