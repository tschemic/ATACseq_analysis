#!/bin/bash

# Preparation and setup of required files

FILES=$(pwd)
WKDIR=$(dirname $FILES)

read -p 'Do you want to retrieve genomic data from the CGD? (yes or no): ' GENEDATA

# Ask for raw data file format

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT

read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW

read -p 'How many threads should be used for the analysis (use 1, if you are not sure): ' THREAD

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
BLKLIST=$FILES/blacklist.bed
mkdir $WKDIR/QC
PICARD=$WKDIR/required_files/picard.jar
mkdir $WKDIR/stats
STATS=$WKDIR/stats

### Fragment sizes for file splitting into nucleosome-free, mono- and dinucleosome fragments
FREESIZE=100
MONOSIZE1=180
MONOSIZE2=240
DISIZE1=315
DISIZE2=437

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
	mkdir $WKDIR/QC_raw
	echo 'Quality control of raw data:'
	if [ $FORMAT == 'bam' ]
	then
		for i in $WKDIR/*.bam
		do
			fastqc -o $WKDIR/QC_raw $i
		done
	else
		for i in $WKDIR/*.fq
		do
			fastqc -o $WKDIR/QC_raw $i
		done
	fi
	multiqc -s -o $WKDIR/QC_raw $WKDIR/QC_raw
else
	echo 'No QC of raw data done.'
fi

# Adapter removal with cutadapt and mapping of all files with NGM

for i in $WKDIR/$(ls $WKDIR | egrep '(\.f.*q$)|(q\.gz$)')
do
	SNAME=$(echo $i | sed 's:/.*/::g')
	cutadapt --interleaved -j $THREAD -q 30 -O 1 -a $ADAPT1 -A $ADAPT2 $i > $i.trimmed.fq.gz 2>$WKDIR/QC/Cutadapt_$SNAME.txt   # removes Illumina TrueSeq adapters from reads (change -a for different adapters); -j specifies number of cores to use, remove if not sure
	rm $i

	ngm -q $i.trimmed.fq.gz -r $GENOME -o $i.trimmed.fq.bam -b -p -Q 30 -t $THREAD # add -p for paired-end data; -t 6 is optional - means 6 threads of the processor are used, if you don't know what to do, remove it; --topn 1 --strata causes ngm to write only uniquely mapping reads to the output
	rm $i.trimmed.fq.gz

	samtools sort -@ $THREAD $i.trimmed.fq.bam -o $i.trimmed.fq.bam.sort.bam   # sort .bam files using samtools
	rm $i.trimmed.fq.bam
	
	# Labelling of duplicated reads and removal of optical duplicates
	java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT I=$i.trimmed.fq.bam.sort.bam O=$i.trimmed.fq.bam.sort.bam.markdup.bam M=$WKDIR/QC/$SNAME.markdup.metrics.txt   ### use REMOVE_SEQUENCING_DUPLICATES=true to remove only optical duplicates
	rm $i.trimmed.fq.bam.sort.bam
	cat $WKDIR/QC/$SNAME.markdup.metrics.txt | sed "s/INPUT=\[.*\]/INPUT=\[$SNAME\]/" > $WKDIR/QC/temp
	mv $WKDIR/QC/temp $WKDIR/QC/$SNAME.markdup.metrics.txt

	bedtools intersect -a $i.trimmed.fq.bam.sort.bam.markdup.bam -b $BLKLIST -v > $i.trimmed.fq.bam.sort.bam.markdup.bam.filt.bam  # removal of reads mapping to mitochondrial loci
	rm $i.trimmed.fq.bam.sort.bam.markdup.bam

	#sambamba view -F "mapping_quality >= 30" -f bam $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam > $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam.qfilt.bam  # not required, already filteres during mapping
	#rm $i.trimmed.fq.bam.sort.bam.filt.bam.markdup.bam
	
	samtools view -H $i.trimmed.fq.bam.sort.bam.markdup.bam.filt.bam > out.sam
	samtools view $i.trimmed.fq.bam.sort.bam.markdup.bam.filt.bam | awk 'BEGIN {FS="\t"; OFS="\t"} {if (($2 == 99)||($2 == 83)||($2 == 147)||($2 == 163)) {print}}' >> out.sam
	samtools view -b -@ $THREAD out.sam > $i.final.bam
	rm $i.trimmed.fq.bam.sort.bam.markdup.bam.filt.bam
	
	samtools sort -n -o $i.final.sortn.bam -@ $THREAD $i.final.bam

	#echo $i >> $WKDIR/QC/flagstat_analysis.txt
	samtools flagstat $i.final.bam >> $WKDIR/QC/$SNAME.txt   # flagstat analysis

	samtools index $i.final.bam

	fastqc -o $WKDIR/QC $i.final.bam
	unzip $WKDIR/QC/$SNAME.final_fastqc.zip
	cat $WKDIR/QC/$SNAME.final_fastqc/fastqc_data.txt | sed 's/Filename	\(.*\).final.bam/Filename	FastQC_\1/' > $WKDIR/QC/$SNAME.final_fastqc/temp
	mv $WKDIR/QC/$SNAME.final_fastqc/temp $WKDIR/QC/$SNAME.final_fastqc/fastqc_data.txt
	zip -um -b $WKDIR/QC/$SNAME.final_fastqc $WKDIR/QC/$SNAME.final_fastqc.zip

done

multiqc -s -o $WKDIR/QC $WKDIR/QC


# Calculation of normalized strand cross-correlation (NSC) and relative strand cross-correlation (RSC)
# CeMM recommendation: NSC should be higher than 1; RSC: the higher the better;
echo "Filename	numReads	estFragLen	corr_estFragLen	PhantomPeak	corr_phantomPeak	argmin_corr	min_corr	Normalized_SCC_(NSC)	Relative_SCC_(RSC)	QualityTag)" > $STATS/spp_results_summary.txt
for i in $WKDIR/*.final.bam
do
Rscript $FILES/run_spp_MT.R -c=$i -savp=$i.plot.pdf -out=$i.sppresults.txt
cat $i.sppresults.txt >> $STATS/spp_results_summary.txt
rm $i.sppresults.txt
mv $i.plot.pdf $STATS
done

# Preparation of coverage files for visualization in IGV

mkdir $WKDIR/IGV_files

for i in $WKDIR/*.final.bam
do
	SNAME=$(echo $i | sed 's:/.*/::g')
	bamCoverage -b $i -o $WKDIR/IGV_files/$SNAME.bw -e -p $THREAD --normalizeUsing CPM
done


# Fragment length splitting and peak calling with nucleosome free fragments
mkdir $WKDIR/macs2_peaks_calling

SQFREESIZE=$(echo "$FREESIZE*$FREESIZE" | bc)
SQMONOSIZE1=$(echo "$MONOSIZE1*$MONOSIZE1" | bc)
SQMONOSIZE2=$(echo "$MONOSIZE2*$MONOSIZE2" | bc)
SQDISIZE1=$(echo "$DISIZE1*$DISIZE1" | bc)
SQDISIZE2=$(echo "$DISIZE2*$DISIZE2" | bc)

for i in $WKDIR/*.final.bam
do
	samtools view -H $i > header.sam
	samtools view $i | awk -v s="$SQFREESIZE" '($9*$9) <= s' > temp1.sam
	samtools view $i | awk -v l="$SQMONOSIZE1" -v u="$SQMONOSIZE2" '($9*$9) > l && ($9*$9) < u' > temp2.sam
	samtools view $i | awk -v l="$SQDISIZE1" -v u="$SQDISIZE2" '($9*$9) > l && ($9*$9) < u' > temp3.sam
	cat header.sam temp1.sam | samtools view -b -@ $THREAD | samtools sort -@ $THREAD > $i.nucfree.bam
	samtools index $i.nucfree.bam
	cat header.sam temp2.sam | samtools view -b -@ $THREAD | samtools sort -@ $THREAD > $i.mononuc.bam
	samtools index $i.mononuc.bam
	cat header.sam temp3.sam | samtools view -b -@ $THREAD | samtools sort -@ $THREAD > $i.dinuc.bam
	samtools index $i.dinuc.bam
	rm *.sam
done

for i in $(ls | grep -E '((nucfree)|(mononuc)|(dinuc)).bam$')
do
	bamCoverage -b $i -o $i.bw -bs 5 -p $THREAD --normalizeUsing CPM
done

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


