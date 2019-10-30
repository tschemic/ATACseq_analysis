### Secondary analysis of ATACseq samples
### Change working directory into the directory with the "final.bam" files
library(tidyverse)
library(RCurl)
library(ggsci)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(ChIPseeker)
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/plot_cleanup.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/convert_to_human_readable_numbers.R", ssl.verifypeer = FALSE)))

### Data import #########################################################################################

gff_file <- "../required_files/C_albicans_SC5314_A22_current_features_haploid_5UTRs.gff"
GeneData <- rtracklayer::import.gff(gff_file)
TxDb.Calbicans.A22.CGD.features <- makeTxDbFromGFF(gff_file, dataSource = "CGD", organism = "Candida albicans SC5314",
                                                   circ_seqs = "Ca22chrM_C_albicans_SC5314")
CGDID_GeneName <- elementMetadata(GeneData)[GeneData$type == "gene", c(5,12)]
CGDID_GeneName$locus_tag <- gsub(pattern = "_", replacement = "", x = CGDID_GeneName$ID)

gff_file2 <- "../required_files/GCF_000182965.3_ASM18296v3_genomic.gff"
GeneData2 <- rtracklayer::import.gff(gff_file2)

CGDID_Entrez <- as.data.frame(elementMetadata(GeneData2[GeneData2$type == "gene", c(20,6)]))
CGDID_Entrez$locus_tag <- gsub(pattern = "CAALFM_", replacement = "", x = CGDID_Entrez$locus_tag)
CGDID_Entrez$Dbxref <- gsub(pattern = "GeneID:", replacement = "", x = CGDID_Entrez$Dbxref)

CaGeneIDs <- base::merge(CGDID_GeneName, CGDID_Entrez, by.x=3, by.y=1, all=TRUE)
colnames(CaGeneIDs) <- c("locus_tag", "CGDID", "GeneName", "EntrezID")

mergedPeaksID <- import("../macs2_peaks_calling/merged_all_peaks.gff")

### Created list of .bam files and sample info ##################################################################
bamFiles <- dir(path = "..", pattern = "final.bam$", full.names = TRUE)
sampleList <- read_tsv("sample_list.txt")
sampleList_wo_contr <- sampleList[sampleList$group != "gDNA", ]
mergedPeaks <- import("../macs2_peaks_calling/merged_all_peaks_annot.bed")
mergedPeaksID <- import("../macs2_peaks_calling/merged_all_peaks.gff")

### Function to plot read numbers mapped to each chromosome ######################################################
plotChromReads <- function(bamList = bamFiles, samplList = sampleList){

  chrom <- tibble(Chromosome = idxstatsBam(bamList[1])[[1]], Length = idxstatsBam(bamList[1])[[2]])
  for (i in 1:dim(samplList)[1]){
    chrom <- cbind(chrom, idxstatsBam(bamList[i])[3])
    colnames(chrom)[i+2] <- samplList[i,2]
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

png(file = "chromReadsPlot.png", res = 600, width = 4200, height = 4200)
plotChromReads()
dev.off()

# Chromosome reads plotted for merged bam files
# Biological replicates (aligned bam files) for ecah condition have been merged using samtools merge
bamFilesMerged <- c("../gDNA_merged_final.bam", "../YPD_merged_final.bam", "../H2O2_merged_final.bam")
sampleListMerged <- tibble(file = c("gDNA_merged_final.bam", "YPD_merged_final.bam", "H2O2_merged_final.bam"), sample = c("gDNA", "YPD", "H2O2"))

chrom_read_plot <- plotChromReads(bamList = bamFilesMerged, samplList = sampleListMerged)
pdf(file = "chromReadsPlot_merged.pdf")
chrom_read_plot
dev.off()
png(filename = "chromReadsPlot_merged.png", res = 600, width = 4200, height = 3600)
chrom_read_plot
dev.off()

### Function to extract the read length distributions of all samples ##############################################

extractFragmentSize <- function(sampList = sampleList){
  insertSizesAll <- read.table(text = "", colClasses = c("numeric", "numeric", "factor"), 
                               col.names = c("insertSizes", "Freq", "sample"))
  for (i in 1:dim(sampList)[1]) {
    fname <- paste0(sampList[[i,1]], ".fragsize")
    sname <- sampList[[i,2]]
    insertSizes <- read_tsv(file = fname, col_names = FALSE)
    insertSizes <- data.frame(table(abs(insertSizes)))
    insertSizes$sample <- sname
    names(insertSizes)[1] <- "insertSizes"
    insertSizes$Freq <- insertSizes$Freq/sum(insertSizes$Freq)
    insertSizesAll <- rbind(insertSizesAll, insertSizes)
  }
  return(insertSizesAll)
}

insertSizesAll <- extractFragmentSize()

### Function to plot the read length distributions of all samples #####################################################

plotFragmentSize <- function(sampList = sampleList[[2]], isizes = insertSizesAll){  
  ### character vector of sample names from sample list can be used as sampList argument to chose specific samples
  
  fragLenPlot <- isizes %>% 
    mutate(InsertSize = as.numeric(as.vector(insertSizes)), Frequency = as.numeric(as.vector(Freq)), 
         Sample = as.factor(as.vector(sample))) %>% 
    .[.$Sample %in% sampList,] %>%
    ggplot(aes(x = InsertSize, y = Frequency, color = Sample)) + 
    geom_line() + cleanup + scale_color_npg() + theme(legend.key = element_blank())
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

# Get fragment sizes for merged bam files and plot them #

fsizes_merged <- extractFragmentSize(bamList = bamFilesMerged, sampList = sampleListMerged)
frag_length_plot_merged <- plotFragmentSize(sampList = sampleListMerged[[2]], isizes = fsizes_merged)

pdf(file = "fragLengthPlot_merged.pdf", height = 5, width = 7.5)
frag_length_plot_merged + scale_color_manual(values = c("grey80", "#e64b35ff", "#4dbbd5ff"))
frag_length_plot_merged + scale_color_manual(values = c("grey80", "#e64b35ff", "#4dbbd5ff")) + 
  scale_y_continuous(trans='log10', limits = c(1e-5, 1e-2)) +
  scale_x_continuous(limits = c(0,800))
dev.off()
png(file = "fragLengthPlot_merged.png", height = 1440, width = 2160, res = 300)
frag_length_plot_merged + scale_color_manual(values = c("grey80", "#e64b35ff", "#4dbbd5ff"))
dev.off()
png(file = "fragLengthPlot_merged_log.png", height = 1440, width = 2160, res = 300)
frag_length_plot_merged + scale_color_manual(values = c("grey80", "#e64b35ff", "#4dbbd5ff")) +
  scale_y_continuous(trans='log10', limits = c(1e-5, 1e-2)) +
  scale_x_continuous(limits = c(0,800))
dev.off()


### Function to calculate number of reads within peaks for each sample ################################################
library(rtracklayer)
library(ChIPQC)

Ca_feat <- import("../required_files/C_albicans_SC5314_A22_current_features_haploid.gff")
CaAnno <- list(version="custom", Ca_features=Ca_feat)
openRegionPeaks <- dir(path = "../macs2_peaks_calling/", pattern = "narrowPeak$", full.names = TRUE)
nucFreeFiles <- dir(path = "..", pattern = "bam.nucfree.bam$", full.names = TRUE)

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

plotRiP <- function(peakCall = peakCallResult) {
  RiPplot <- peakCall %>% mutate(RiP. = RiP.*100) %>% 
    gather(key = "Data", value = "Peaks", RiP.:Peaks) %>%
    ggplot() + geom_col(aes(x = Sample, y = Peaks, fill = Data), position = "dodge") +
    scale_y_continuous(sec.axis = sec_axis(~./100, name = "% of reads in peaks"), name = "No. of peaks") +
    labs(x = "") + coord_flip() + cleanup + scale_fill_npg(name = "", labels = c("No. of peaks", "% of reads in peaks"))
}

RiPplot <- plotRiP()
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

annotatePeaks <- function(peaks = mergedPeaks){
  anno <- annotatePeak(peaks, TxDb = TxDb.Calbicans.A22.CGD.features, tssRegion = c(-2000, 0), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
  return(anno)
}

mergedPeaksAnno <- annotatePeaks()

plotAnnoPie(mergedPeaksAnno)
plotAnnoBar(mergedPeaksAnno)
upsetplot(mergedPeaksAnno)
vennpie(mergedPeaksAnno)
#upsetplot(mergedPeaksAnno, vennpie=TRUE) # not working
plotDistToTSS(mergedPeaksAnno, title = "Distribution of peaks relative to TSS", 
              ylab = "Peaks (%) (5'->3')")
covplot(mergedPeaks)
#covplot(test_peaks, weightCol = "score")  ### peaks over all chromosomes - including weightcol works only with single sample - requires some value to plot (e.g. score column in MACS2 files)

#promoter <- getPromoters(TxDb=TxDb.Calbicans.A22.CGD.features, upstream=500, downstream=200)
#tagMatrix <- getTagMatrix(mergedPeaks, windows=promoter)
#plotAvgProf(tagMatrix, xlim=c(-750, 200), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
#tagHeatmap(tagMatrix, xlim=c(-500, 200), color="red") ### does not work propely

#plotAvgProf2(mergedPeaks, TxDb=TxDb.Calbicans.A22.CGD.features, upstream=1000, downstream=200,
 #            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

pdf(file = "Peak_annotation_plots.pdf", height = 5)
upsetplot(mergedPeaksAnno)
plotAnnoBar(mergedPeaksAnno)
plotAnnoPie(mergedPeaksAnno)
dev.off()

png(file = "upsetplot_peak_annotation.png", res = 300, width = 1920, height = 1440)
upsetplot(mergedPeaksAnno)
dev.off()
png(file = "barplot_peak_annotation.png", res = 300, width = 1920, height = 800)
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



### Differential ATACseq analysis with edgeR and annotation of peaks to genes ####################################
library(ChIPseeker)
library(edgeR)

pval = 0.05 # set pvalue threshold (only relevant for plotting)
pval_adjust = "BH"  # set method for p-value adjustment for multiple testing (= FDR calculation)
cutoff = c(-1,1)  # set log2 fold change line in plots (e.g. "c(-1,1)" means lines will be drawn at 2-fold up and down)

annotMergedPeaks <- annotatePeak(mergedPeaksID, TxDb = TxDb.Calbicans.A22.CGD.features, tssRegion = c(-2000, 0))

#targets <- readTargets(file = "Targets.txt") # this imports the data in the Targets.txt file = experiment setup

targets <- as.data.frame(sampleList_wo_contr[,c(1,3,2)])
targets[1] <- paste0("../", targets[[1]], ".nucfree.bam.sortn.bam.count.txt.crop.txt")
targets[2] <- as.factor(targets[[2]])
colnames(targets) <- c("files", "group", "description")
#levels(targets$condition) <- c(control, treatment)


d <- readDGE(targets, header = FALSE) # this reads in the count data
d$genes <- data.frame(peakID = row.names(d$counts))
m <- match(d$genes$peakID, annotMergedPeaks@anno$ID)
d$genes$geneID <- annotMergedPeaks@anno$geneId[m]
d$genes$annotation <- annotMergedPeaks@anno$annotation[m]
m <- match(d$genes$geneID, CaGeneIDs$CGDID)
d$genes$GeneName <- CaGeneIDs$GeneName[m]
row.names(d$samples) <- d$samples$description
d$genes$EntrezID <- CaGeneIDs$EntrezID[m]

keep <- filterByExpr(d) ### filtering low count peaks - from current edgeR manual
d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis
plotMDS(d, labels = d$samples$description, col=as.numeric(d$samples$group))

png(filename = "MDSplot.png", res = 300, width = 1920, height = 1440)
plotMDS(d, labels = d$samples$description, col=as.numeric(d$samples$group))
dev.off()

# Model fitting without batch effect correction ##########################################################################
# create design matrix
design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)
my.contrasts <- makeContrasts(
  H2O2vsYPD = H2O2-YPD,
  levels=design)

# Estimating the dispersions and plot them
d <- estimateDisp(d, design, robust=TRUE) ## does both dispersion estimations in one step - suggested in edgeR manual; add design matrix for multifactor experiments

pdf("disp.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
plotBCV(d)
dev.off()# closes and saves the pdf file with the plot

fit <- glmQLFit(d, design, robust = TRUE)
plotQLDisp(fit)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,c("H2O2vsYPD")])
#qlf <- glmQLFTest(fit, contrast=c(1,-1,0,0,0,0)) ### same as above
topTags(qlf)
summary(decideTests(qlf))

plotMD(qlf)
abline(h=c(-1, 1), col="blue")

res <- as.data.frame(topTags(qlf, n=Inf))


# Model fitting with batch effect correction ############################################################################
Cond <- factor(targets$group)
Cond <- relevel(Cond, ref = "YPD")
Repl <- factor(gsub(pattern = "(YPD)|(H2O2)", replacement = "repl", targets$description))

design_paired <- model.matrix(~Repl+Cond)
rownames(design_paired) <- colnames(d)
#logFC <- predFC(d,design,prior.count=1,dispersion=0.05)
#cor(logFC[,4:6])

dp <- estimateDisp(d, design_paired, robust=TRUE)
pdf(file = "disp.pdf")
plotBCV(dp)
dev.off()

fitp <- glmQLFit(dp, design_paired, robust = TRUE)
plotQLDisp(fitp)

#qlfp <- glmQLFTest(fitp, coef = "Strainhir1:Timet8") # alternative for testing for one condition
qlfp <- glmQLFTest(fitp, coef = "CondH2O2")
#qlfp <- glmQLFTest(fitp, contrast = c(0,0,0,0,0,0,0,-1,1)) # to test for time point specificity
topTags(qlfp)
summary(decideTests(qlfp))
plotMD(qlfp)

res_paired <- as.data.frame(topTags(qlfp, n=Inf))
rng <- as.data.frame(ranges(annotMergedPeaks@anno))
info <- elementMetadata(annotMergedPeaks@anno)[c(5,14,8,9)]
sqn <- as.data.frame(annotMergedPeaks@anno@seqnames)
rngInfo <- cbind(rng,info,sqn)
res_coord <- merge(res_paired, rngInfo, by.x=1, by.y=4, all.x=TRUE)
bed <- as.data.frame(res_coord[res_coord$FDR < 0.05 , 
                               #(res_coord$logFC > 0), 
                               c(17,15,16,2)])
write_tsv(unique(bed), "Sign_up_genes.bed", col_names = FALSE)
#sign_prom_bed <- as.data.frame(res_coord[res_coord$annotation == "Promoter (<=1kb)" & res_coord$width < 1000 & res_coord$logFC > 0 & res_coord$FDR < 0.05, c(17,11,12,1,6,2)])
#strnd <- GeneData[GeneData$type == "gene",] %>%  as.data.frame()
#strnd <- strnd[,c(10,5)]
#sign_prom_bed <- merge(sign_prom_bed, strnd, by.x=6, by.y=1, all.x=TRUE) # export of peak regions for downstream analysis e.g. motif search
#write_tsv(sign_prom_bed[,c(2:7)], "Sign_up_peaks_t4_1000bp_P1kb.bed", col_names = FALSE)




############################################################################################################

d <- readDGE(targets, header = FALSE) # this reads in the count data

# Filter low expression tags (cpm<1)
keep <- rowSums(cpm(d)> 1) >= 3 # this creates an R object with the genes/transcripts that have more than 1 counts per million - change the 1 if you want to modify the filter; the 3 at the end is the group size - change this if you have more or less replicates than 3
d <- d[keep,] # this filters out all genes/transcripts with less than 1 counts per million
d$samples$lib.size <- colSums(d$counts) # this calculates the library size again after the filtering step

# Normalization (TMM)
d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis

# Estimating the dispersions and plot them
d <- estimateCommonDisp(d, verbose=FALSE) # this calculates the common dispersion over all genes and samples - is used during the diff. expr. anaylsis
d <- estimateTagwiseDisp(d) # this calculates the dispersion for each gene/transcript - is used during the diff. expr. anaylsis

pdf("disp.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
plotBCV(d)
dev.off()# closes and saves the pdf file with the plot

# Differential expression, lines with # signs are just control commands to check the analysis

et <- exactTest(d, pair=c('YPD','H2O2')) # this command performs the differential expression analysis

#detags <- rownames(topTags(et, n=20))
result <- summary(de <- decideTestsDGE(et, p=pval, adjust=pval_adjust)) # this command gives the number of differentially expressed genes; add "lfc = 1" for logFC -1/1 as additional cutoff

detags <- rownames(d)[as.logical(de)] # this command saves the row names of differentially expressed genes as "detags"

# generate smear plots (cloud of dots in log2CPM~log2FC axes)

pdf('Diff_scatterplot.pdf')# change the name of the .pdf file; this opens a pdf file to store the next plot
plotSmear(et, de.tags=detags, xlab = "Average log2 counts per million", ylab = "log2 fold change") # this plots the analysis results; use "panel.first = NULL" to remove grid
abline(h = cutoff, col = "blue")
dev.off() # closes and saves the pdf file with the plot 

# Export results
res <- as.data.frame(topTags(et, n=Inf)) # this saves the fold change table in "tt"

annotMergedPeaks <- annotatePeak(mergedPeaksID, TxDb = TxDb.Calbicans.A22.CGD.features, tssRegion = c(-3000,0))
annotMergedPeaksDF <- as.data.frame(annotMergedPeaks)
res <- merge(res, annotMergedPeaksDF, by.x=0, by.y=10, all.x=TRUE)
res <- merge(res, CaGeneIDs, by.x=21, by.y=2, all.x=TRUE)

### Create Vulcano plot ########################
pdf(file = "Vulcano_plot_diff_edgeR.pdf")
ggplot(res, aes(logFC, -log10(FDR))) + geom_point(alpha=0.4) + cleanup + 
geom_hline(yintercept = -log10(pval), colour='red') + geom_vline(xintercept = cutoff, colour='blue') + 
xlab('Fold change (log2)') + ylab('Adjusted p-value (-log10)') +
geom_text(aes(label = ifelse(res$FDR < 0.05 & 
        (res$logFC < -1 | res$logFC > 1.1), res$GeneName, "")))
dev.off()


####### GO term enrichment analysis ##############################################################

library(AnnotationHub)
library(clusterProfiler)

CalbicansEntrez <- AnnotationHub()[['AH73663']]

resProm <- res[res$annotation == "Promoter (<=1kb)",c(1,2,3,6,23,25,26)]
resPromMax <- resProm %>% filter(!is.na(resProm$EntrezID) & resProm$FDR < 0.05) %>% dplyr::select(logFC, EntrezID) %>% 
  group_by(EntrezID) %>% 
  summarise(max = ifelse(abs(max(logFC)) > abs(min(logFC)), max(logFC), min(logFC)))

#upgenes <- resProm[resProm$FDR < 0.05 & resProm$logFC > 0,]
upgenes <- resPromMax[resPromMax$max > 0,]
#downgenes <- resProm[resProm$FDR < 0.05 & resProm$logFC < 0,]
downgenes <- resPromMax[resPromMax$max < 0,]
#geneList <- unique(na.omit(upgenes$EntrezID))
geneList <- upgenes$EntrezID
foldChange <- as.vector(upgenes$max)
names(foldChange) <- upgenes$EntrezID

go <- enrichGO(geneList, OrgDb = CalbicansEntrez, ont = "BP", readable = TRUE, qvalueCutoff = 0.05)
go_simpl <- simplify(go, cutoff = 0.7)
dplot <- dotplot(go_simpl, showCategory =20) + 
  scale_color_gradient(name = "Adj. p-value", low = "#cb181d", high = "#fee0d2") +
  labs(size = "No. of Genes")
dplot

cplot <- cnetplot(go_simpl, showCategory = 15,  foldChange = foldChange) +
  scale_color_gradient(low = "#ccece6", high = "#00441b") +
  labs(color = "Fold Change\n(log2)", size = "No. of Genes")
cplot

pdf(file = "GOdplot_diff_up.pdf", height = 6, width = 9)
dplot
dev.off()
png(filename = "GOdplot_diff_up.png", res = 150, width = 1440, height = 920)
dplot
dev.off()

pdf(file = "GOcnetplot_diff_up.pdf", width = 10)
cplot
dev.off()
png(filename = "GOcnetplot_diff_up.png", res = 150, width = 1920, height = 920)
cplot
dev.off()

### Gene set enrichment analysis using GO terms as gene sets ###########################################
gseData <- resProm %>% filter(!is.na(resProm$EntrezID)) %>% dplyr::select(logFC, EntrezID) %>% 
  group_by(EntrezID) %>% 
  summarise(max = ifelse(abs(max(logFC)) > abs(min(logFC)), max(logFC), min(logFC)))
gseData <- gseData[order(gseData$max, decreasing = TRUE),]
gseList <- gseData$max
names(gseList) <- gseData$EntrezID

gsego <- gseGO(gseList, ont = "BP", OrgDb = CalbicansEntrez, pvalueCutoff = 0.05, exponent = 0.5)
gsego_simpl <- simplify(gsego, cutoff = 0.5)
dotplot(gsego_simpl, showCategory = 20)
cnetplot(gsego_simpl)


### Overlay ATAC-seq with RNAseq data (H2O2 treated log phase cells in YPD) #############################
### adjust accordingly for other datasets to overlay

RNAseq <- read_tsv("RNAseqResults_NGMedgeRcpm1_wt.csv")

A22orf19IDs <- read_tsv("A22_orf19_ID_table.txt", col_names = FALSE)
names(A22orf19IDs) <- c("CGDID", "orf19ID")

RNAseq <- merge(RNAseq, A22orf19IDs, by.x=1, by.y=2)

atacRNAseq <- merge(res, RNAseq, by.x=1, by.y=5)

RNAseqPlot <- ggplot(data = atacRNAseq, aes(x=logFC, y=-log10(FDR), color=`wtu-wtt_logFC`)) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype="dashed") +
  geom_point() +
  scale_color_gradient2(low = "steelblue1", mid = "white", high = "red") +
  labs(x="ATACseq fold change (log2)", y="ATACseq adj. p-value\n(-log10)", color="RNAseq fold change\n(log2)") +
  cleanup

pdf("RNAseq_compare_vulcano.pdf", height = 4)
RNAseqPlot
dev.off()

png("RNAseq_compare_vulcano.png", res = 300, height = 920, width = 1920)
RNAseqPlot
dev.off()

### Export RNAseq regulated genes with coordinates ###
#genes <- as.data.frame(GeneData[elementMetadata(GeneData)$type == "gene",]) # old version
genes <- as.data.frame(genes(TxDb.Calbicans.A22.CGD.features))
RNAseq_coord <- merge(RNAseq, genes, by.x=5, by.y=6)
up4x <- RNAseq_coord[RNAseq_coord$`wtu-wtt_logFC` > 2 & RNAseq_coord$`wtu-wtt_FDR` < 0.05, c(6,7,8,1)]
down4x <- RNAseq_coord[RNAseq_coord$`wtu-wtt_logFC` < -2 & RNAseq_coord$`wtu-wtt_FDR` < 0.05, c(6,7,8,1)]
write_tsv(up4x, "../coverage_analysis/RNAseq_up4x.bed", col_names = FALSE)
write_tsv(down4x, "../coverage_analysis/RNAseq_down4x.bed", col_names = FALSE)


### Plot ATACseq peaks with fold change H2O2 vs YPD across chromosomes ########################################################
### template from: https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneExpression/GeneExpression.html ###
library(karyoploteR)

resPromMaxCGDID <- resProm %>% dplyr::select(geneId, logFC, GeneName, Row.names) %>%
  group_by(geneId, GeneName) %>% 
  summarise(max = ifelse(abs(max(logFC)) > abs(min(logFC)), max(logFC), min(logFC))) %>% as.data.frame()

ca.genes <- genes(TxDb.Calbicans.A22.CGD.features)
meta <- merge(elementMetadata(ca.genes), resPromMaxCGDID, by.x=1, by.y=1, all.x=TRUE)
row.names(meta) <- meta$gene_id
mcols(ca.genes) <- meta[names(ca.genes), c("gene_id", "GeneName", "max")]
ordered <- ca.genes[order(ca.genes$max, na.last = TRUE, decreasing = TRUE),]
CaChrom <- import("../required_files/chrMA.bed")
kp <- plotKaryotype(genome = CaChrom[seqnames(CaChrom) != "Ca22chrM_C_albicans_SC5314"])
kp <- kpPlotMarkers(kp, ordered[1:20], labels = elementMetadata(ordered)$GeneName[1:20], text.orientation = "horizontal")

filtered.ca.genes <- ca.genes[!is.na(ca.genes$max)]
fc.ymax <- ceiling(max(abs(range(filtered.ca.genes$max))))
fc.ymin <- -fc.ymax
kp <- plotKaryotype(genome = CaChrom[seqnames(CaChrom) != "Ca22chrM_C_albicans_SC5314"])
kp <- kpPoints(kp, data = filtered.ca.genes, y=filtered.ca.genes$max, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
#kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin)


### Create coverage plot with locally scaled signal from matrix obtained from deeptolls computeMatrix #################
# deepTools used to get matrix of coverage values at ATGs with merged bigWig files (all replicates merged):
# computeMatrix reference-point -S YPD_nucfree_merged.bw YPD_mononuc_merged.bw -R CaA22_all_genes.bed -o YPD_nucfree_mononuc_all_genes.npz --outFileNameMatrix YPD_nucfree_mononuc_all_genes.npz.tab -a 200 -b 1000 -bs 5 -p 5

# TSS data obtained from: http://www.yeastss.org/jbrowse/JBrowse_data/Candida_albicans/CalbconsensusClusters.txt
# NCBI chromosome names have been replaced and the modified file can be found in "/required_files"

TSSs <- read_tsv("../required_files/CalbconsensusClusters_CGDID.txt", col_names = FALSE) # import TSS (transcription start site data)
names(TSSs) <- c("clusterID", "chrom", "start", "stop", "strand", "dominant_ctss", "tpm", "tpm.dominant_ctss", "q_0.1", "q_0.9", "interquantile_width") # for more info on dataset see www.yeastss.org
TSSs_filt <- TSSs[TSSs$tpm.dominant_ctss > 5, c(2,6,5)]
TSSs_filt$dominant_ctss_stop <- ifelse(TSSs_filt$strand == "+", TSSs_filt$dominant_ctss + 1, TSSs_filt$dominant_ctss)
TSSs_filt[TSSs_filt$strand == "-",]$dominant_ctss <- TSSs_filt[TSSs_filt$strand == "-",]$dominant_ctss - 1
TSSs_filt$name <- "."
TSSs_filt$score <- "0"

write.table(TSSs_filt[,c(1,2,4,5,6,3)], "../coverage_analysis/CalbconsensusClusters_CGDID_filt_gr5.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep ="\t")


### Import of coverage matrix files produced by deepTools ###################################
coverage_data_file <- "../coverage_analysis/YPD_H2O2_nucfree_mononuc_filt_TSSs_a1000.npz.tab"  # coverage data from deepTools
upstream_bp <- -1000  # bp upstream of TSS used in deepTools
downstream_bp <- 1000  # bp downstream of TSS used in deepTools
bs <- 5  # bin size used in deepTools
csamples <- c("Ynf", "Ymn", "Hnf", "Hmn")  # short sample names of coverage samples (bigWig files) used for deepTools
csample_names <- c("YPD Nucleosome-free", "YPD Nucleosome-occupied", "H2O2 Nucleosome-free", "H2O2 Nucleosome-occupied")  # full sample names of coverage samples (bigWig files) used for deepTools


dt_matrix <- read_tsv(coverage_data_file, skip = 3, col_names = FALSE)
xvalues <- seq(upstream_bp, downstream_bp - bs, bs)
bin_number <- dim(dt_matrix)[2] / length(csamples)

cnames <- NULL
for (i in 1:length(csamples)) {
  n <- paste0(csamples[i], "_", seq(1, bin_number))
  cnames <- c(cnames, n)
}

colnames(dt_matrix) <- cnames

coveragePlotData <- function(x = dt_matrix, y = csamples){
  cov_data <- read.table(text = "", colClasses = c("numeric", "numeric", "character"), 
                               col.names = c("coord", "coverage", "sample"))
  for(i in 1:length(y)) {
    sampl_name <- y[i] # gets the sample name
    dt_matrix_sampl <- dplyr::select(x, starts_with(sampl_name)) # selects one sample
    dt_matrix_sampl[is.na(dt_matrix_sampl)] <- 0 # assigns 0 to missing values
    dt_matrix_max <- apply(dt_matrix_sampl, 1, max) # calculates the maximum values for each row (=locus)
    dt_matrix_sampl <- dt_matrix_sampl[dt_matrix_max != 0,] # removes rows (=loci) with no signal across the whole region
    dt_matrix_sampl_sc <- apply(dt_matrix_sampl, 1, function(x) x/max(x)) # calculates maximum for each row (=locus)
    cov_means <- (rowMeans(dt_matrix_sampl_sc) - min(rowMeans(dt_matrix_sampl_sc))) / (max(rowMeans(dt_matrix_sampl_sc)) - min(rowMeans(dt_matrix_sampl_sc)))
    dt_matrix_sampl_mean <- data.frame(coord = xvalues, coverage = cov_means, sample = sampl_name) # calculates the mean for each bin
    cov_data <- rbind(cov_data, dt_matrix_sampl_mean) # merges data for all samples into 1 data frame
  }
  return(cov_data)
}

cov_plot_data <- coveragePlotData()

### Plotting of coverage data ###########
cov_plot <- ggplot(cov_plot_data, aes(x = coord, y = coverage, color = sample)) +
  geom_line(lwd = 1) + cleanup + scale_color_brewer(labels = csample_names, palette = "Dark2") +
  labs(x = "Postion relative to TSS (bp)", y = "Fraction of signal", color = "") +
  theme(legend.key = element_blank(), legend.position = c(0.15, 0.9))

cov_plot_anno <- cov_plot + geom_vline(aes(xintercept = 0), lty=2, color = "grey") +
  geom_segment(aes(x = 0, y = 1.03, xend = 100, yend = 1.03), arrow = arrow(length = unit(0.2, "cm")), color = "grey") +
  geom_text(aes(x=55, y=1.04, label = "TSS"), color = "grey", vjust = -0.4, fontface = "plain", size = 3)
cov_plot_anno

### File export #############  
pdf(file = "../coverage_analysis/CovPlot_filt_exactpos_TSS_a1000_RNAseq4xup.pdf", height = 4)
cov_plot_anno
dev.off()

png(filename = "../coverage_analysis/CovPlot_filt_exactpos_TSS_a1000_RNAseq4xup.png", res = 300, height = 960, width = 1920)
cov_plot_anno
dev.off()

############ Coverage plot from deepTools coverage matrix without local scaling ######################
### import of coverage data is described above ###
coveragePlotDataRaw <- function(x = dt_matrix, y = csamples){
  cov_data <- read.table(text = "", colClasses = c("numeric", "numeric", "character"), 
                         col.names = c("coord", "coverage", "sample"))
  for(i in 1:length(y)) {
    sampl_name <- y[i] # gets the sample name
    dt_matrix_sampl <- dplyr::select(x, starts_with(sampl_name)) # selects one sample
    dt_matrix_sampl <- dt_matrix_sampl[complete.cases(dt_matrix_sampl),]
    cov_means <- colMeans(dt_matrix_sampl)
    dt_matrix_sampl_mean <- data.frame(coord = xvalues, coverage = cov_means, sample = sampl_name) # calculates the mean for each bin
    cov_data <- rbind(cov_data, dt_matrix_sampl_mean) # merges data for all samples into 1 data frame
  }
  return(cov_data)
}

cov_plot_data_raw <- coveragePlotDataRaw()

### Plotting of coverage data ###########
cov_plot <- ggplot(cov_plot_data_raw, aes(x = coord, y = coverage, color = sample)) +
  geom_line(lwd = 1) + cleanup + scale_color_brewer(labels = csample_names, palette = "Dark2") +
  labs(x = "Postion relative to TSS (bp)", y = "Relative read coverage (log2)", color = "") +
  theme(legend.key = element_blank(), legend.position = c(0.15, 0.9))

cov_plot_anno <- cov_plot + geom_vline(aes(xintercept = 0), lty=2, color = "grey") +
  geom_segment(aes(x = 0, y = max(cov_plot_data_raw$coverage)*1.2, xend = 100, yend = max(cov_plot_data_raw$coverage)*1.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "grey") +
  geom_text(aes(x=55, y=max(cov_plot_data_raw$coverage)*1.25, label = "TSS"), color = "grey", vjust = -0.4, fontface = "plain", size = 3)
cov_plot_anno

### File export #############  
pdf(file = "../coverage_analysis/CovPlot_.pdf", height = 4)
cov_plot_anno
dev.off()

png(filename = "../coverage_analysis/CovPlot_.png", res = 300, height = 960, width = 1920)
cov_plot_anno
dev.off()

#################### Subsetting TSS data for H2O2-regulated genes #########################

RNAseq <- read_tsv("RNAseqResults_NGMedgeRcpm1_wt.csv")
A22orf19IDs <- read_tsv("A22_orf19_ID_table.txt", col_names = FALSE)
names(A22orf19IDs) <- c("CGDID", "orf19ID")
RNAseq <- merge(RNAseq, A22orf19IDs, by.x=1, by.y=2)
RNAseq_filt <- RNAseq[RNAseq$`wtu-wtt_logFC` > 2 & RNAseq$`wtu-wtt_FDR` < 0.05,]
RNAseq_filt2 <- RNAseq[RNAseq$`wtu-wtt_logFC` < 0.58 & RNAseq$`wtu-wtt_logFC` > -0.58 & RNAseq$`wtu-wtt_FDR` > 0.05, ]

TSSs_filt_GR <- import.bed("../coverage_analysis/CalbconsensusClusters_CGDID_filt_gr5.bed")
annotatePeaks <- function(peaks = merged_peaks){
  anno <- annotatePeak(peaks, TxDb = TxDb.Calbicans.A22.CGD.features, sameStrand = TRUE)
  return(anno)
}

  #TSSs_anno <- annotatePeaks(peaks = TSSs_filt_GR)
  #TSSs_anno_DF <- as.data.frame(TSSs_anno@anno[TSSs_anno@anno$annotation == "Promoter (<=1kb)",])

  #TSSs_RNAseq <- merge(TSSs_anno_DF, RNAseq_filt, by.x=14, by.y=5)
  #TSSs_RNAseq_bed <- TSSs_RNAseq[,c(2,3,4,1,8,6)]
  #TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "+",]$end <- TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "+",]$end + 1
  #TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "-",]$start <- TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "-",]$start - 1


mergeTSSs <- function(tss = TSSs_filt_GR, rnaseq = RNAseq_filt){
  TSSs_anno <- annotatePeaks(peaks = tss)
  TSSs_anno_DF <- as.data.frame(TSSs_anno@anno[TSSs_anno@anno$annotation == "Promoter (<=1kb)",])
  TSSs_RNAseq <- merge(TSSs_anno_DF, rnaseq, by.x=14, by.y=5)
  TSSs_RNAseq_bed <- TSSs_RNAseq[,c(2,3,4,1,8,6)]
  TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "+",]$end <- TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "+",]$end + 1
  TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "-",]$start <- TSSs_RNAseq_bed[TSSs_RNAseq_bed$strand == "-",]$start - 1
  return(TSSs_RNAseq_bed)
}

tss_RNAseq_4xup <- mergeTSSs()

tss_RNAseq_unchanged <- mergeTSSs(rnaseq = RNAseq_filt2)


### Exported file can be used in deepTools to calculate coverage at H2O2-induced genes
write.table(tss_RNAseq_4xup, "../coverage_analysis/CalbconsensusClusters_CGDID_filt_gr5_RNAseq4xup.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep ="\t")
write.table(tss_RNAseq_unchanged, "../coverage_analysis/CalbconsensusClusters_CGDID_filt_gr5_RNAsequnchanged.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep ="\t")


### Cut site visualization ####################################################################
library(MotifDb)
library(Biostrings)
library(GenomicAlignments)
library(soGGi)

binding_sites <- "motif_search/Cap1_binding_sites.bed"

nucFreeFiles <- dir(path = "..", pattern = "bam.nucfree.bam$", full.names = TRUE)
fimo_hits <- import(binding_sites)

extractCutSites <- function(reads){
  read1 <- first(reads)
  read2 <- second(reads)
  Firsts <- resize(granges(read1), fix = "start", 1)
  First_Pos_toCut <- shift(granges(Firsts[strand(read1) == "+"]), 4) # extracts only forward reads (only forward read is used because one insertion event produces two reads, 1 forward and one reverse)
  #First_Neg_toCut <- shift(granges(Firsts[strand(read1) == "-"]), -5)

  Seconds <- resize(granges(read2), fix = "start", 1)
  Second_Pos_toCut <- shift(granges(Seconds[strand(read2) == "+"]), 4)
  #Second_Neg_toCut <- shift(granges(Seconds[strand(read2) == "-"]), -5)

  test_toCut <- c(First_Pos_toCut, Second_Pos_toCut)
  cutsCoverage <- coverage(test_toCut)
  return(cutsCoverage)
}

covPlotData <- function(bamFiles = nucFreeFiles, rangs = fimo_hits, flanks = 500){  
  for (i in 1:length(bamFiles)) {
    GApair <- readGAlignmentPairs(bamFiles[i])
    cutsCov <- extractCutSites(GApair)
    sname <- sampleList_wo_contr[[2]][i]
    cap1CutsSample <- regionPlot(cutsCov, testRanges = rangs, style = "point", samplename = sname,
                         format = "rlelist", distanceAround = flanks)
    if (i == 1) {
      Cap1_cuts_list <- list(cap1CutsSample)
    } else {
      Cap1_cuts_list <- c(Cap1_cuts_list, list(cap1CutsSample))
    }
  }
  Cap1_cuts <- c(Cap1_cuts_list[[1]], Cap1_cuts_list[[2]], Cap1_cuts_list[[3]], Cap1_cuts_list[[4]], Cap1_cuts_list[[5]], Cap1_cuts_list[[6]]) # adjust this to number of samples
  return(Cap1_cuts)
}


Cap1_cuts_open <- covPlotData()

regPlot <- plotRegion(Cap1_cuts_open, outliers = 0.001, colourBy = "Sample", groupBy="Sample") + 
  ggtitle("NucFree Cuts Centred on Cap1") +
  labs(x = "") +
  theme(legend.key = element_rect(fill = NA)) +
  scale_color_npg() + cleanup
regPlot

pdf(file = "RegionPlot.pdf", height = 5, width = 7)
regPlot
dev.off()

png(filename = "RegionPlot.png", res = 300, height = 1440, width = 2200)
regPlot
dev.off()


###################### Import and GO term enrichment with deepTools coverage data after clustering #######
cl_regions_file <- "../coverage_analysis/log2ratio_H2O2_YPD_merged_transcripts_k4.npz.png.bed"
cl_regions <- import.bed(cl_regions_file)
gene_anno <- annotatePeak(cl_regions, TxDb = TxDb.Calbicans.A22.CGD.features, tssRegion = c(-100, 0), sameStrand = TRUE)
gene_annoDF <- as.data.frame(gene_anno)
gene_annoDF <- gene_annoDF[,c(1,16,17,20,18,5,6,9)]
gene_annoDF_names <- merge(gene_annoDF, CaGeneIDs, by.x=4, by.y=2, all.x=TRUE, sort = FALSE)

# merging with H2O2 RNAseq data
RNAseq <- read_tsv("../downstream_analysis/RNAseqResults_NGMedgeRcpm1_wt.csv")
A22orf19IDs <- read_tsv("../downstream_analysis/A22_orf19_ID_table.txt", col_names = FALSE)
names(A22orf19IDs) <- c("CGDID", "orf19ID")
RNAseq <- merge(RNAseq, A22orf19IDs, by.x=1, by.y=2)
cluster_RNAseq_merged <- merge(gene_annoDF_names, RNAseq, by.x=1, by.y=5, all.x=TRUE)

clusterDF <- cluster_RNAseq_merged[cluster_RNAseq_merged$NA. == "cluster_1" &
                                     !is.na(cluster_RNAseq_merged$EntrezID) &
                                     !is.na(cluster_RNAseq_merged$`wtu-wtt_logFC`),]
cluster <- clusterDF$EntrezID
foldChange <- clusterDF$`wtu-wtt_logFC`
names(foldChange) <- clusterDF$EntrezID

go <- enrichGO(cluster, OrgDb = CalbicansEntrez, ont = "BP", readable = TRUE, qvalueCutoff = 0.05)
go_simpl <- simplify(go, cutoff = 0.7)
dplot <- dotplot(go, showCategory =20) + 
  scale_color_gradient(name = "Adj. p-value", low = "#cb181d", high = "#fee0d2") +
  labs(size = "No. of Genes")
dplot

cplot <- cnetplot(go, showCategory = 15,  foldChange = foldChange) +
  scale_color_gradient2(low = "steelblue1", mid = "grey95", high = "red") +
  labs(color = "RNAseq fold Change\n(log2)", size = "No. of Genes")
cplot

pdf(file = "../coverage_analysis/GOdplot_cluster1.pdf", height = 6, width = 9)
dplot
dev.off()
pdf(file = "../coverage_analysis/GOcplot_cluster1.pdf", height = 9, width = 15)
cplot
dev.off()

png(filename = "../coverage_analysis/GOdplot_cluster1.png", res = 300, width = 2880, height = 1840)
dplot
dev.off()
png(filename = "../coverage_analysis/GOcplot_cluster1.png", res = 200, width = 2880, height = 1840)
cplot
dev.off()


