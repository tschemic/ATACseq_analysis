#!/bin/bash

# Preparation and setup of required files

WKDIR=$(grep "working directory" ./config_file.txt | cut -d ":" -f 2)
FILES=$WKDIR/required_files

read -p 'Do you want to retrieve genomic data from the CGD? (yes or no): ' GENEDATA

# Ask for raw data file format

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT

read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW

# read -p 'Are the libraries prepared in a strand-specific way? (yes or no): ' STRANDED   ### only for RNAseq data

if [ $GENEDATA == 'yes' ]
then
	echo 'Retrieving genomic data from CGD.'

	wget http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz  ## include WKDIR/required_files folder in wget!!!
	wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff
	cat C_albicans_SC5314_A22_current_features.gff | egrep -v "Ca22chr[1-7R]B" > C_albicans_SC5314_A22_current_features_haploid.gff
	gunzip C_albicans_SC5314_A22_current_chromosomes.fasta.gz
	cat C_albicans_SC5314_A22_current_chromosomes.fasta | egrep ">Ca22chr[1-7RM][A_]" | sed 's/>//g' | sed 's/(/1	/g' | sed 's/ nucleotides)//g' | sed 's/ /	/g' > chrMA.bed
	bedtools getfasta -fi C_albicans_SC5314_A22_current_chromosomes.fasta -bed chrMA.bed | fold -w 60 | sed 's/:1-[0-9]*//g' > C_albicans_SC5314_A22_current_chromosomesAM.fasta
else
	echo 'No genomic data are retrieved.'
fi


GENOME=$WKDIR/required_files/C_albicans_SC5314_A22_current_chromosomesAM.fasta
FEATURES=$WKDIR/required_files/C_albicans_SC5314_A22_current_features_haploid.gff
ADAPT1=$(cat $WKDIR/required_files/config_file.txt | grep Read1: | cut -d ":" -f 2)
ADAPT2=$(cat $WKDIR/required_files/config_file.txt | grep Read2: | cut -d ":" -f 2)
# rRNA=$WKDIR/required_files/Ca_A22chrAM_rRNAloci.bed  ### only for RNAseq data
MITO=$FILES/chrM.bed
mkdir $WKDIR/QC
PICARD=$(cat $WKDIR/required_files/config_file.txt | grep "picard file path:" | cut -d ":" -f 2)
mkdir $WKDIR/stats
STATS=$WKDIR/stats

if [ $FORMAT == 'bam' ]
then
	echo 'File format is bam.'
	for i in $WKDIR/*.bam
	do
		bamToFastq -i $i -fq $i.fq
	done
elif [ $FORMAT == 'fastq' ]
then
	echo 'File format is fastq.'
else
	echo 'Invalid file format! Options are "bam" or "fastq".'
	exit
fi

# QC of raw data


if [ $QCRAW == 'yes' ]
then
	echo 'Quality control of raw data:'
	if [ $FORMAT == 'bam' ]
	then
		for i in $WKDIR/*.bam
		do
			fastqc -o $WKDIR/QC $i
		done
	else
		for i in $WKDIR/*.fq
		do
			fastqc -o $WKDIR/QC $i
		done
	fi
else
	echo 'No QC of raw data done.'
fi

# Adapter removal with cutadapt and mapping of all files with NGM

for i in $WKDIR/*.fq
do
	SNAME=$(echo $i | sed 's:/.*/::g')
	cutadapt --interleaved -j 5 -q 30 -a $ADAPT1 -A $ADAPT2 $i > $i.trimmed.fq.gz 2>$WKDIR/QC/$SNAME.cutadapt.report.txt   # removes Illumina TrueSeq adapters from reads (change -a for different adapters); -j specifies number of cores to use, remove if not sure
	#rm $i

	ngm -q $i.trimmed.fq.gz -r $GENOME -o $i.trimmed.fq.bam -b -t 5 -p --topn 1 --strata # add -p for paired-end data; -t 6 is optional - means 6 threads of the processor are used, if you don't know what to do, remove it; --topn 1 --strata causes ngm to write only uniquely mapping reads to the output
	#rm $i.trimmed.fq.gz

	samtools sort $i.trimmed.fq.bam -o $i.trimmed.fq.bam.sort.bam   # sort .bam files using samtools
	#rm $i.trimmed.fq.bam

	bedtools intersect -a $i.trimmed.fq.bam.sort.bam -b $MITO -v > $i.trimmed.fq.bam.sort.bam.filt.bam  # removal of reads mapping to mitochondrial loci
	#rm $i.trimmed.fq.bam.sort.bam

	#samtools rmdup -s $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam.rmdup.bam  # removal of duplicated reads

	# Labelling of duplicated reads and removal of optical duplicates
	java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT I=$i.trimmed.fq.bam.sort.bam.filt.bam O=$i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam M=$WKDIR/QC/$SNAME.markdup.metrics.txt   ### use REMOVE_SEQUENCING_DUPLICATES=true to remove only optical duplicates
	#rm $i.trimmed.fq.bam.sort.bam.filt.bam

	sambamba view -F "mapping_quality >= 30" -f bam $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam > $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam.qfilt.bam
	#rm $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam

	#echo $i >> $WKDIR/QC/flagstat_analysis.txt
	samtools flagstat $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam.qfilt.bam >> $WKDIR/QC/$SNAME.flagstat_analysis.txt   # flagstat analysis

	samtools index $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam.qfilt.bam

	fastqc -o $WKDIR/QC $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam.qfilt.bam

done

multiqc -s -o $WKDIR/QC $WKDIR/QC


# Calculation of normalized strand cross-correlation (NSC) and relative strand cross-correlation (RSC)
# CeMM recommendation: NSC should be higher than 1; RSC: the higher the better;
echo "Filename	numReads	estFragLen	corr_estFragLen	PhantomPeak	corr_phantomPeak	argmin_corr	min_corr	Normalized_SCC_(NSC)	Relative_SCC_(RSC)	QualityTag)" > $STATS/spp_results_summary.txt
for i in $WKDIR/*.qfilt.bam
do
Rscript $FILES/run_spp_MT.R -c=$i -savp=$i.plot.pdf -out=$i.sppresults.txt
cat $i.sppresults.txt >> $STATS/spp_results_summary.txt
rm $i.sppresults.txt
mv $i.plot.pdf $STATS
done

# Preparation of coverage files for visualization in IGV

mkdir $WKDIR/IGV_files

for i in $WKDIR/*.qfilt.bam
do
	SNAME=$(echo $i | sed 's:/.*/::g')
	bamCoverage -b $i -o $WKDIR/IGV_files/$SNAME.bw -e -p 5 --normalizeUsing CPM
done


# peak calling
# CeMM default settings: --nomodel --extsize 147
#mkdir $WKDIR/macs2_peaks_calling
#for i in $WKDIR/*.qfilt.bam
#do
#SNAME=$(echo $i | sed 's:/.*/::g' | cut -d "." -f 1)  ### this extracts the sample name from the bam file name - adjust accordingly, otherwise replace "$SNAME" in script with "$i" (without quotes)
#macs2 callpeak --nomodel --extsize 147 -g 14521502 -n $SNAME -t $i --outdir $WKDIR/macs2_peaks_calling  ### 14521502 bp is the size of the haploid A22 C. albicans genome without mitochondria
#done

# Converts .narrowPeak file from peak calling to .gff files
#for i in $WKDIR/macs2_peaks_calling/*.narrowPeak
#do
#cp $i $WKDIR/input.txt
#$WKDIR/narrowPeak_gff_conversion.r
#echo "##gff-version 3" > $i.gff
#cat $WKDIR/output.txt >> $i.gff
#rm $WKDIR/input.txt
#rm $WKDIR/output.txt
#done

# counts number of reads in called peaks
#for i in $WKDIR/*.rmdup.bam
#do
#SNAME=$(echo $i | sed 's:/.*/::g' | cut -d "." -f 1)
#htseq-count -f bam -s no -t peak -i ID $i $WKDIR/macs2_peaks_calling/$SNAME"_peaks.narrowPeak.gff" > $WKDIR/macs2_peaks_calling/$SNAME.count.txt
#done

# Calculates the number of reads within peaks
#for i in $WKDIR/macs2_peaks_calling/*.count.txt
#do
#echo $i >> $WKDIR/macs2_peaks_calling/FRiP.txt
#echo "Reads in peaks:" >> $WKDIR/macs2_peaks_calling/FRiP.txt
#cat $i | grep -v ^_ | awk '{sum += $2} END {print sum}' >> $WKDIR/macs2_peaks_calling/FRiP.txt
#PEAKS=$(cat $i | grep -v ^_ | awk '{sum += $2} END {print sum}')
#echo "Total reads:" >> $WKDIR/macs2_peaks_calling/FRiP.txt
#cat $i | awk '{sum += $2} END {print sum}' >> $WKDIR/macs2_peaks_calling/FRiP.txt
#TOTAL=$(cat $i | awk '{sum += $2} END {print sum}')
#echo "FRiP:" >> $WKDIR/macs2_peaks_calling/FRiP.txt
#echo $PEAKS/$TOTAL | bc -l >> $WKDIR/macs2_peaks_calling/FRiP.txt
#done

#####################################################################################################################################################################################



# Required files for analysis (these files must be in the directory specified below as "WKDIR"):
# 1) raw read data files either in .bam format or fastq format (files must have .bam or .fq extension)
# 2) haploid C. albicans genome (e.g. A chromosomes of assembly 22 - mitochondrial chromosome can be included) in .fasta format (specifiy filename in line 11, if necessary))
# 3) fasta file with adaptor sequences for trimming (e.g. nextera_wo_duplicates.fa; specifiy filename in line 12, if necessary)
# 4) .bed file with ID, start and stop coordinates of mitochondrial chromosome (e.g. chrM.bed; specifiy filename in line 13, if necessary)
# 5) Chromosome length file (e.g. Ca_chromosomesA.bed) containing chromosome lengths
# 6) Script files: analysis script (ATACseq_script.sh), run_spp_MT.R and narrowPeak_gff_conversion.r files

#WKDIR=/home/tschemic/Data/ATACseq/190208_ATACseq_pipeline_construction
#GENOME=$WKDIR/C_albicans_SC5314_version_A22-s07-m01-r73_chromosomesA_haploid.fasta
#ADAPT=$WKDIR/nextera_wo_duplicates.fa
#MITO=$WKDIR/chrM.bed
#CHROM=$WKDIR/Ca_chromosomesA.bed
#mkdir $WKDIR/stats
#STATS=$WKDIR/stats

# FastQC quality control of libraries
#for i in $WKDIR/*.bam
#do
#mkdir $WKDIR/fastqc
#fastqc $i -o $WKDIR/fastqc
#done

#multiqc $WKDIR/fastqc

# conversion from .bam to .fastq format (required by skewer - next step)
#for i in $WKDIR/*.bam
#do
#bedtools bamtofastq -i $i -fq $i.fq

# trimming of adapters - skewer replaced by cutadapt in script version 2.1
#skewer -x $ADAPT $i.fq
#rm $i.fq
#mv $i'-trimmed.log' $STATS

#cutadapt -a CTGTCTCTTATACACATCT -o $i'-trimmed.fq' $i.fq

# Mapping of all .bam files with NGM
#ngm -q $i'-trimmed.fastq' -r $GENOME -o $i.trimmed.fastq.sam --top 1 --strata # -t 6 is optional - means 6 threads of the processor are used, if you don't know what to do, remove it
											# use "--top 1 --strata" to allow only uniquely mapped reads (in practice (e.g. reads with mapping quality 60) NextGenMap only finds one good alignment anyway
#rm $i'-trimmed.fastq'

# convert .sam to .bam using samtools (only required if mapping output is .sam)
#samtools view -bS $i.trimmed.fastq.sam | samtools sort > $i.trimmed.fastq.sam.sorted.bam
#rm $i.trimmed.fastq.sam

# filter out mitochondrial reads
#bedtools intersect -v -a $i.trimmed.fastq.sam.sorted.bam -b $MITO > $i.trimmed.fastq.sam.sorted.bam.mitominus.bam

# Calculate stats for
#SNAME=$(echo $i | sed 's:/.*/::g' | cut -d "." -f 1)  ### this extracts the sample name from the bam file name - adjust accordingly, otherwise replace "$SNAME" in script with "$i" (without quotes)

#echo "Reads in $i:" > $STATS/$SNAME.stats
#samtools view $i | wc -l >> $STATS/$SNAME.stats
#echo "Mapped reads in $i.trimmed.fastq.sam.sorted.bam:" >> $STATS/$SNAME.stats
#samtools flagstat $i.trimmed.fastq.sam.sorted.bam >> $STATS/$SNAME.trimmed.fastq.sam.sorted.bam.flagstat
#cat $STATS/$SNAME.trimmed.fastq.sam.sorted.bam.flagstat | grep "mapped (" >> $STATS/$SNAME.stats
#echo "Mapped reads in $i.trimmed.fastq.sam.sorted.bam.mitominus.bam:" >> $STATS/$SNAME.stats
#samtools flagstat $i.trimmed.fastq.sam.sorted.bam.mitominus.bam | grep "mapped (" >> $STATS/$SNAME.stats
#rm $i.trimmed.fastq.sam.sorted.bam

# filter aligned files for mapping quality higher than 30 (phred score)
#sambamba view -F "mapping_quality >= 30" -f bam $i.trimmed.fastq.sam.sorted.bam.mitominus.bam > $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam
#echo "Mapped reads in $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam:" >> $STATS/$SNAME.stats
#samtools flagstat $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam | grep "mapped (" >> $STATS/$SNAME.stats
#rm $i.trimmed.fastq.sam.sorted.bam.mitominus.bam

# index bam files
#samtools index $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam

# remove duplicated reads
#samtools rmdup -s $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam.rmdup.bam
#echo "Mapped reads in $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam.rmdup.bam:" >> $STATS/$SNAME.stats
#samtools flagstat $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam.rmdup.bam | grep "mapped (" >> $STATS/$SNAME.stats

# index bam files
#samtools index $i.trimmed.fastq.sam.sorted.bam.mitominus.bam.filtered.bam.rmdup.bam

#done

## Calculation of normalized strand cross-correlation (NSC) and relative strand cross-correlation (RSC)
## CeMM recommendation: NSC should be higher than 1; RSC: the higher the better;
#echo "Filename	numReads	estFragLen	corr_estFragLen	PhantomPeak	corr_phantomPeak	argmin_corr	min_corr	Normalized_SCC_(NSC)	Relative_SCC_(RSC)	QualityTag)" > $STATS/spp_results_summary.txt
#for i in $WKDIR/*.rmdup.bam
#do
#Rscript $WKDIR/run_spp_MT.R -c=$i -savp=$i.plot.pdf -out=$i.sppresults.txt
#cat $i.sppresults.txt >> $STATS/spp_results_summary.txt
#rm $i.sppresults.txt
#mv $i.plot.pdf $STATS
#done



## Creates bigWig files for visualization e.g. in IGV
#mkdir $WKDIR/IGV_files
#for i in $WKDIR/*.filtered.bam
#do
#SNAME=$(echo $i | sed 's:/.*/::g')
#genomeCoverageBed -ibam $i -bg -g $CHROM > $WKDIR/IGV_files/$SNAME.fullread.bedGraph
#READS=$(samtools view $i | wc -l)
#awk -v reads="${READS}" '{{ $4 = ($4 / reads); print}}' $WKDIR/IGV_files/$SNAME.fullread.bedGraph > $WKDIR/IGV_files/$SNAME.fullread.bedGraph.norm
#genomeCoverageBed -ibam $i -5 -bg -g $CHROM > $WKDIR/IGV_files/$SNAME.bedGraph
#awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' $WKDIR/IGV_files/$SNAME.bedGraph $WKDIR/IGV_files/$SNAME.bedGraph > $WKDIR/IGV_files/$SNAME.bedGraph.norm #normalization factor=1000000, change if desired; it will only change the absolute numbers for the coverage track in IGV
#rm $WKDIR/IGV_files/$SNAME.fullread.bedGraph
#rm $WKDIR/IGV_files/$SNAME.bedGraph
#done

#for i in $WKDIR/*.rmdup.bam
#do
#SNAME=$(echo $i | sed 's:/.*/::g')
#genomeCoverageBed -ibam $i -bg -g $CHROM > $WKDIR/IGV_files/$SNAME.fullread.bedGraph
#READS=$(samtools view $i | wc -l)
#awk -v reads="${READS}" '{{ $4 = ($4 / reads) * 100000; print}}' $WKDIR/IGV_files/$SNAME.fullread.bedGraph > $WKDIR/IGV_files/$SNAME.fullread.bedGraph.norm #normalization factor=100000, change if desired; it will only change the absolute numbers for the coverage track in IGV
#genomeCoverageBed -ibam $i -5 -bg -g $CHROM > $WKDIR/IGV_files/$SNAME.bedGraph
#awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' $WKDIR/IGV_files/$SNAME.bedGraph $WKDIR/IGV_files/$SNAME.bedGraph > $WKDIR/IGV_files/$SNAME.bedGraph.norm #normalization factor=1000000, change if desired; it will only change the absolute numbers for the coverage track in IGV
#rm $WKDIR/IGV_files/$SNAME.fullread.bedGraph
#rm $WKDIR/IGV_files/$SNAME.bedGraph
#done

#for i in $WKDIR/IGV_files/*.norm
#do
#bedGraphToBigWig $i $CHROM $i.bw
#rm $i
#done
