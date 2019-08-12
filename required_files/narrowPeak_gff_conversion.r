#!/usr/bin/Rscript

rdata <- read.table("merged_all_peaks_annot.bed",
                    sep="\t", quote = "")
rdata <- rdata[,c(1,2,3)]
rdata$source <- "macs2"
rdata$type <- "peak"
rdata$score <- "."
rdata$strand <- "."
rdata$phase <- "."
#rdata$temp1 <- "ID=peak_"
rdata$attributes <- paste("ID=peak_", seq(1, length(rdata$source), 1), sep="")
rdata <- rdata[,c(1,4,5,2,3,6,7,8,9)]
write.table(rdata, file="merged_all_peaks.gff", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
