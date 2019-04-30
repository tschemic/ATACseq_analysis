# ATACseq_analysis

Script for analysis of ATAC-seq data obtained from Candida albicans. This script is used in the Kuchler lab (http://cdl.univie.ac.at/) at MFPL (https://www.mfpl.ac.at/de.html).
It is based on this ATACseq pipeline: https://github.com/epigen/open_pipelines/blob/master/pipelines/atacseq.md

This repository conatains a pipeline for the primary analysis of Illumina short read sequencing ATACseq data (paired-end) obtained from the fungal pathogen Candida albicans. It includes the retrieval of genomic data required for analysis from the Candida Genome Database (CGD, http://www.candidagenome.org/), quality control of the raw data, trimming and mapping of reads, removal of reads mapping to mitochonria and removing duplicates. Only the ATACseq raw data in .bam or .fastq format (compressed or uncompressed) have to be provided by the user.

#Tools required for analysis:

samtools (http://www.htslib.org/)

bedtools (https://bedtools.readthedocs.io/en/latest/)

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

MultiQC (https://multiqc.info/)

cutadapt (https://cutadapt.readthedocs.io/en/stable/)

NextGenMap (https://github.com/Cibiv/NextGenMap/wiki)

sambamba (http://lomereiter.github.io/sambamba/)

DeepTools (https://deeptools.readthedocs.io/en/develop/)

R spp package (https://cran.r-project.org/web/packages/spp/index.html)

All the above-mentioned tools have to be included in yout PATH environment.

Picard tools (https://broadinstitute.github.io/picard/) Download the picard.jar file and add the path to that file to the config_file.txt

#Usage:

Clone the repository and copy the raw data into the ATACseq_analysis directory.

Adjust the path to the ATACseq_analysis directory in the config_file.txt.

Adjust the path to your picard.jar file in the config_file.txt.

Change the adapter sequence for read trimming in the config_file.txt if necessary. By default it contains the Illumina Nextera adapter. For adapter sequences see: https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-10.pdf 
Change into the required_files directory and run the analysis script (analysis_script.sh).
