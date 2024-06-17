#!/usr/bin/bash

#-----README-----#
# Mike Martinez 2023, Rosenberg Lab, UCHC
#----------------#

# The purpose of this script is to merge fastq files from two lanes of a flow cell (sample replicates)
# Run this script if you have 2 forward and 2 reverse reads for a single sample (F/R pair from one flow cell lane and a F/R pair from another flow cell lane)
# This script will generate a directory within the directory that contains your fastq files called `mergedReplicates`
# Within this directory, the forward reads from both lanes will be merged into one file, and likewise for the reverse reads. 
#----------------#

#-----Let's start by listing all of the files present in the directory
allSamples=(*.fastq.gz)

echo -e "Listing all samples present in directory\n"
for sample in "${allSamples[@]}"; do
    echo ${sample}
done

#----Now we want to isolate all the reads that came from one sequencing run, in this case, we have reads from lane L003 and L004

# Initialize two empty arrays
FirstLane=()
FirstLaneSamples=()
SecondLane=()
SecondLaneSamples=()

# Iterate through the samples and append to respective array
for sample in "${allSamples[@]}"; do
    if [[ $sample == *L003* ]]; then
        FirstLane+=(${sample})
    elif [[ $sample == *L004* ]]; then
        SecondLane+=(${sample})
    fi
done

# Obtain the base sample name for each sample
for sample in "${FirstLane[@]}"; do 
    basename=$(echo "${sample}" | cut -d'_' -f1)
    FirstLaneSamples+=(${basename})
done
echo -e "\n"

for sample in "${SecondLane[@]}"; do
    basename=$(echo "${sample}" | cut -d'_' -f1)
    SecondLaneSamples+=(${basename})
done
echo -e "\n"

# Check to see if the number of samples is the same between the two arrays. 
if [ ${#FirstLane[@]} -eq ${#SecondLane[@]} ]; then
    echo "Arrays are of equal length."
else
    echo "Error! Two arrays are not of equal length"
    exit 1
fi

# Check to see if the elements are the same between the two arrays. 
if [ "${FirstLaneSamples[*]}" == "${SecondLaneSamples[*]}" ]; then
    echo "The arrays have the same elements."
else
    echo "Error! Elements are not the same between the two arrays."
    exit 2
fi

#-----Concatenate all the forward and reverse reads

# Initialize 4 arrays to hold the forward and reverse reads from the first lane and second lane respectively
FL_forward=()
FL_forwardNames=()
FL_reverse=()
Fl_reverseNames=()
SL_forward=()
SL_forwardNames=()
SL_reverse=()
SL_reverseNames=()

echo "Segregating forward and reverse reads from Lane 1"
for sample in "${FirstLane[@]}"; do
    if [[ $sample == *R1* ]]; then
        FL_forward+=(${sample})
        basename=$(echo "${sample}" | cut -d'_' -f1)
        FL_forwardNames+=(${basename})
    elif [[ $sample == *R2 ]]; then
        FL_reverse+=(${sample})
        basename=$(echo "${sample}" | cut -d'_' -f1)
        FL_reverseNames+=(${basename})
    fi
done

echo "Segregating forward and reverse reads from Lane 2"
for sample in "${SecondLane[@]}"; do
    if [[ $sample == *R1* ]]; then
        SL_forward+=(${sample})
        basename=$(echo "${sample}" | cut -d'_' -f1)
        SL_forwardNames+=(${basename})
    elif [[ $sample == *R2* ]]; then
        SL_reverse+=(${sample})
        basename=$(echo "${sample}" | cut -d'_' -f1)
        SL_reverseNames+=(${basename})
    fi
done

#----- Concatenate all the forward reads and reverse reads from lane 1 and lane 2 into one sample


# Create a directory to hold the merged replicates
mkdir -p ./mergedReplicates

# Iteratively merge the replicates and rename
for ((i=0; i<${#FL_forward[@]}; i++)); do
    cat ${FL_forward[i]} ${SL_forward[i]} >> ./mergedReplicates/${FL_forwardNames[i]}_R1_Merged.fastq.gz
    cat ${FL_reverse[i]} ${SL_reverse[i]} >> ./mergedReplicates/${SL_reverseNames[i]}_R2_Merged.fastq.gz
done






