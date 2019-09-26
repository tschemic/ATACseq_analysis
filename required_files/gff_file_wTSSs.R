library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(GenomicFeatures)

# import TSS data
tsss <- read_tsv("CalbconsensusClusters_CGDID.txt", col_names = FALSE)
names(tsss) <- c("clusterID",	"chr",	"start",	"end",	"strand",	"dominant_ctss",	"score",	"score.dominant_ctss",	"q_0.1",	"q_0.9",	"interquantile_width")

# Select TSSs with highest signal for each gene
tsss_bed <- tsss[,c(2,6,8,5)]
tsss_bed$end <- tsss_bed$dominant_ctss
tsss_bed[tsss_bed$strand == "+",]$end <- tsss_bed[tsss_bed$strand == "+",]$end + 1
tsss_bed[tsss_bed$strand == "-",]$dominant_ctss <- tsss_bed[tsss_bed$strand == "-",]$dominant_ctss - 1
tsss_bed <- as.data.frame(tsss_bed[,c(1,2,5,4,3)])
names(tsss_bed) <- c("chrom", "start", "end", "strand", "score")
tsss_GR <- makeGRangesFromDataFrame(tsss_bed, keep.extra.columns = TRUE)

gff_file = "C_albicans_SC5314_A22_current_features_haploid.gff"
annotatePeaks <- function(peaks = tsss_GR, gff = gff_file){
  TxDb.Calbicans.A22.CGD.features <- makeTxDbFromGFF(gff, dataSource = "CGD", organism = "Candida albicans",
                                                     circ_seqs = "Ca22chrM_C_albicans_SC5314")
  anno <- annotatePeak(peaks, TxDb = TxDb.Calbicans.A22.CGD.features, sameStrand = TRUE, tssRegion = c(-3000, 0))
  return(anno)
}

tsss_anno <- annotatePeaks(tsss_GR)
tsss_anno_DF <- as.data.frame(tsss_anno)
tsss_anno_DF <- tsss_anno_DF[grepl(x = tsss_anno_DF$annotation, pattern = "Promoter"),]

genes_max_score <- tsss_anno_DF %>% group_by(geneId) %>% summarise(chrom = seqnames[which.max(score)],
                                                                   start = start[which.max(score)], 
                                                                   end = end[which.max(score)],
                                                                   score = max(score), strand = strand[which.max(score)])
genes_max_score <- genes_max_score[,c(2,3,4,1,5,6)]
names(genes_max_score)[4] <- "name"
#write.table(genes_max_score, "../coverage_analysis/CalbconsensusClusters_CGDID_max_score.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep ="\t")

# import gff file
gff <- read_tsv(gff_file, col_names = FALSE, skip = 9)
names(gff) <- c("chrom","source","type","start","end","score","strand","phase","attributes")

# modify gff file to include 5' UTRs
exons1 <- gff[grepl(x = gff$attributes, pattern = "-T-E1"),]
exons1 <- exons1 %>% separate(col = attributes, into = c("geneID", "rest"), sep = "-T-E1;", remove = FALSE)
exons1$geneID <- gsub(pattern = "ID=", replacement = "", x = exons1$geneID)

exons1_tss <- merge(exons1, genes_max_score, by.x=10, by.y=4)
exons1_tss[exons1_tss$strand.x == "+",]$start.x <- exons1_tss[exons1_tss$strand.x == "+",]$start.y
exons1_tss[exons1_tss$strand.x == "-",]$end.x <- exons1_tss[exons1_tss$strand.x == "-",]$end.y
exons1_tss <- exons1_tss[,c(2,3,4,5,6,7,8,9,10)]
names(exons1_tss) <- c("chrom","source","type","start","end","score","strand","phase","attributes")

non_exons1 <- gff[!grepl(x = gff$attributes, pattern = "-T-E1"),]
exons1_non_tss <- merge(exons1, genes_max_score, by.x=10, by.y=4, all=TRUE)
exons1_non_tss <- exons1_non_tss[is.na(exons1_non_tss$chrom.y),]
exons1_non_tss <- exons1_non_tss[,c(2,3,4,5,6,7,8,9,10)]
names(exons1_non_tss) <- c("chrom","source","type","start","end","score","strand","phase","attributes")

RNA <- non_exons1[grepl(pattern = "RNA", x = non_exons1$type),]
RNA <- RNA %>% separate(col = attributes, into = c("geneID", "rest"), sep = "-T;", remove = FALSE)
RNA$geneID <- gsub(pattern = "ID=", replacement = "", x = RNA$geneID)
RNA_tss <- merge(RNA, genes_max_score, by.x=10, by.y=4)
RNA_tss[RNA_tss$strand.x == "+",]$start.x <- RNA_tss[RNA_tss$strand.x == "+",]$start.y
RNA_tss[RNA_tss$strand.x == "-",]$end.x <- RNA_tss[RNA_tss$strand.x == "-",]$end.y
RNA_tss <- RNA_tss[,c(2,3,4,5,6,7,8,9,10)]
names(RNA_tss) <- c("chrom","source","type","start","end","score","strand","phase","attributes")
RNA_non_tss <- merge(RNA, genes_max_score, by.x=10, by.y=4, all=TRUE)
RNA_non_tss <- RNA_non_tss[is.na(RNA_non_tss$chrom.y),]
RNA_non_tss <- RNA_non_tss[,c(2,3,4,5,6,7,8,9,10)]
names(RNA_non_tss) <- c("chrom","source","type","start","end","score","strand","phase","attributes")
non_exons1_nonRNA <- non_exons1[!grepl(pattern = "RNA", x = non_exons1$type),]

final_gff <- rbind(non_exons1_nonRNA, RNA_tss, RNA_non_tss, exons1_tss, exons1_non_tss)
final_gff <- final_gff[order(final_gff$chrom, final_gff$start),]

# export new gff file; add gff3 header in linux command line or text editor
write_tsv(final_gff, "C_albicans_SC5314_A22_current_features_haploid_5UTRs_woHeader.gff", col_names = FALSE)
