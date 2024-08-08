#!/bin/bash
#SBATCH --job-name=RNA_Mapping
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=80G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --array=[0-9]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=micmartinez@uchc.edu
#SBATCH --output RNA_Mapping-%j.out
#SBATCH --error RNA_Mapping-%j.err



echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

#Make a temporary directory
# ADJUST FILE PATHS AS NEEDED
mkdir -p /labs/Rosenberg/mmartinez/temp
export TMPDIR=/labs/Rosenberg/mmartinez/pipelineTest

## Current GenomeVersion : Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
##           GTF         :  Rattus_norvegicus.mRatBN7.2.105.gtf


#Specify the paths to the reference genome files and adapters
# ADJUST PATHS AS NEEDED
ALIGNERINDEX="/labs/Rosenberg/mmartinez/HSap/Homo_sapiens"
GTF="/labs/Rosenberg/mmartinez/Homo_sapiens.GRCh38.105.gtf"
SPLICE_SITE="/labs/Rosenberg/mmartinez/splice_site"
ADAPTERFILE="/isg/shared/apps/Trimmomatic/0.39/adapters/NexteraPE-PE.fa"

#Path where the raw fastq files live
# ADJUST PATHS AS NEEDED
#DataDir=/labs/Rosenberg/Shan/Walnut_Study_RNAseq/data/raw_fastq
DataDir=/labs/Rosenberg/Shan/Walnut_Study_RNAseq/data/raw_fastq_test_run/

#Create a merged-reads directory
	#The -p flag ensures that the directory does not already exist
mkdir -p 01_merged_reads

#Move into the newly created directory
cd ${DataDir}

#Extract file basename
SAMPLES=($(ls -1 *L00*_R1*.gz | cut -d"_" -f1 ))

#Access the current sample and $SLURM_ARRAY_TASK_ID specifies the current job task
sampleID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

#Find the forward and reverse reads associated with a sample name and store them in lexicographical order
R1files=$(ls ${sampleID}*R1* | sort -)
R2files=$(ls ${sampleID}*R2* | sort -)

echo -e "\n Order and ID of R1 Files : ${R1files} \n\n Order and ID R2 Files : ${R2files} \n\n"

#Concatenate the forward reads
cat ${R1files} >> /labs/Rosenberg/mmartinez/01_merged_reads/${sampleID}_R1.fastq.gz
#Concatenate the reverse reads
cat ${R2files} >> /labs/Rosenberg/mmartinez/01_merged_reads/${sampleID}_R2.fastq.gz

#Back up two directories. Need to explicitly outline the file path here!!!!!
cd /labs/Rosenberg/mmartinez/
#Create a new directory for QC
mkdir ./02_qualityQC

##########################################################

#Move into the QC directory
cd ./02_qualityQC

# Create two subdirectories in the 02_qualityQC directory to hold paired and unpaired reads respectively. 
mkdir -p ./TrimReads
mkdir -p ./singles

# Load the Trimmomatic module
module load Trimmomatic/0.39

# Run Trimmomatic
java -jar $Trimmomatic PE -threads 6 \
        ../01_merged_reads/${sampleID}_R1.fastq.gz \
        ../01_merged_reads/${sampleID}_R2.fastq.gz \
        ./TrimReads/trim_${sampleID}_R1.fastq.gz ./singles/${sampleID}_singles.fastq.gz \
        ./TrimReads/trim_${sampleID}_R2.fastq.gz ./singles/${sampleID}_singles.fastq.gz \
        ILLUMINACLIP:${ADAPTERFILE}:2:30:10:5:true \
        SLIDINGWINDOW:4:25 MINLEN:45 HEADCROP:2

# Move into the TrimReads (paired) directory
cd ./TrimReads

#Check if the trimmed file exists and is not emptyclear.  If it is not empty, remove the merged fastq file to clear space.
if [[ -f trim_${sampleID}_R1.fastq.gz && -s trim_${sampleID}_R1.fastq.gz ]]; then
	echo "TRIMMOMATIC : ${sampleID} Trim file exists and not empty"
	rm ../../01_merged_reads/${sampleID}_R*.fastq.gz
else
	echo "ERROR: TRIMMOMATIC: ${sampleID} trimmed files do not exist or are empty"
        exit 1
fi

# Back up to the 02_QualityQC directory
cd ..

# Drop the Trimmomatic module and load the fastqc module
module unload Trimmomatic/0.39
module load fastqc

# Make a subdirectory in the 02_qualityQC directory
mkdir -p TRIMfastqc_OUT

# Run fastqc specifying the new TRIMfastqc_OUT directory as the output location. These are html files so they won't take up too much space. 
fastqc -t 6 -o ./TRIMfastqc_OUT ./TrimReads/trim_${sampleID}_R1.fastq.gz ./TrimReads/trim_${sampleID}_R2.fastq.gz

# Do the same for the raw unpaired fastq files. 
mkdir -p RAWfastqc_OUT

fastqc -t 6 -o ./RAWfastqc_OUT ../01_merged_reads/${sampleID}_R1.fastq.gz ../01_merged_reads/${sampleID}_R2.fastq.gz

#####
module purge

# In the main directory, create a new directory for mapping and move into it.
mkdir -p ../03_mapping
cd ../03_mapping

# If a bam file already exists, do nothing. If not, start mapping.
if [ -e "${sampleID}_mapped_sort.bam" ]; then
        echo "${sampleID}_mapped_sort.bam   file EXISTS"
else
        module load hisat2/2.2.1

        mkdir ${sampleID}_tmp

        hisat2 -p 6 --known-splicesite-infile ${SPLICE_SITE} \
        -x ${ALIGNERINDEX} \
        -1 ../02_qualityQC/TrimReads/trim_${sampleID}_R1.fastq.gz \
        -2 ../02_qualityQC/TrimReads/trim_${sampleID}_R2.fastq.gz \
        -S ./${sampleID}_tmp/${sampleID}.sam
        
        # Now that we ran the mapping, we can remove the trimmed fastq files to save space. 
        rm /labs/Rosenberg/mmartinez/02_qualityQC/TrimReads/trim_${sampleID}_R1.fastq.gz
		rm /labs/Rosenberg/mmartinez/02_qualityQC/TrimReads/trim_${sampleID}_R2.fastq.gz

        cd ./${sampleID}_tmp

        module load samtools/1.9
        
        # Convert the SAM file to a binary alignment map
        samtools view -@ 6 -bhS ${sampleID}.sam -o ${sampleID}_mapped.bam
        
        # If the binary alignment file exists and is not empty/clear, remove the original SAM file and sort the BAM file
        if [[ -f ${sampleID}_mapped.bam && -s ${sampleID}_mapped.bam ]]; then
                echo "${sampleID}_mapped_sort.bam exist and not empty"
                rm ${sampleID}.sam
				samtools sort -@ 6 ${sampleID}_mapped.bam -o ${sampleID}_mapped_sort.bam
				
		# If the mapped and sorted BAM file exists and it not empty/clear, remove the mapped unsorted BAM file to save space.
		if [[ -f ${sampleID}_mapped_sort.bam && -s ${sampleID}_mapped_sort.bam ]]; then
                	echo "${sampleID}_mapped_sort.bam exist and not empty"
                	rm ${sampleID}_mapped.bam
		else
			echo "${sampleID}_mapped_sort.bam not exist or empty"
			exit 1
		fi



        else
                echo "${sampleID}_mapped_sort.bam not exist or empty";
                exit 1
        fi

		# Get some stats on the mapped sorted BAM file
        mv ${sampleID}_mapped_sort.bam ../
        samtools flagstat ${sampleID}_mapped_sort.bam
        cd ../

fi


module purge

mkdir -p ../04_counts

cd ../04_counts

module load htseq/0.11.0

htseq-count -s reverse -r pos -t exon -i gene_id -f bam ../03_mapping/${sampleID}_mapped_sort.bam ${GTF} > ${sampleID}.counts

if [[ -f ${sampleID}.count && -s ${sampleID}.count ]]; then
	echo "${sampleID}.counts exist and not empty"
	rm ../03_mapping/${sampleID}_tmp/${sampleID}_mapped_sort.bam


module purge


date
