#!/bin/bash

paired="/home/user/Desktop/RNAseq_analysis/sample_files/fastq_files/trimmo/pair_read"
unpaired="/home/user/Desktop/RNAseq_analysis/sample_files/fastq_files/trimmo/unpair_read"

index_path="/home/user/Desktop/RNAseq_analysis/sample_files/BWA_index/genome_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
aligned_file="/home/user/Desktop/RNAseq_analysis/sample_files/alignmed_Output"

file=$1
steps_count=0
while read -r line; do
    echo "Processing the SRR number  $line "
    echo -e "\n"
     
 #Prefetch     
    echo  "************************************** Prefetch process started**********************************************"
    echo "prefetch -v -p  $line  -O ~/Desktop/RNAseq_analysis/sample_files/prefetch_files/$line" 
    prefetch -p -v  $line -o ~/Desktop/RNAseq_analysis/sample_files/prefetch_files/$line.sra
    echo -e "\n "
    echo -e"/n Prefetch process for $line has  completed"
    echo -e "\n "
  
  #fasterq dump 
     echo  "************************************** fasterq-dump  process started**********************************************"
     
     echo "fasterq-dump ~/Desktop/RNAseq_analysis/sample_files/prefetch_files/$line.sra -p -O ~/Desktop/RNAseq_analysis/sample_files/fastq_files -t ~/Desktop/RNAseq_analysis/sample_files/temp -e 6 "
    
    fasterq-dump ~/Desktop/RNAseq_analysis/sample_files/prefetch_files/$line.sra -p -O ~/Desktop/RNAseq_analysis/sample_files/fastq_files -t ~/Desktop/RNAseq_analysis/sample_files/temp -e 6  
    echo -e "\n "
    echo -e"/n fastewrq-dump  process has  completed"
    echo -e "\n "
    echo " "
    echo "**************************************  $line process completed **********************************************"
    echo -e "\n "
    echo -e "\n "
    echo -e "\n "
    
    
    # Trimmometric process   

    echo "Trimmometric process begin"
    echo -e "\n"
    echo  "************************************** Trimmometry process started for $line *********************************************"
    echo -e "\n" 
    echo "trimmomatic PE -threads 6   ~/Desktop/RNAseq_analysis/sample_files/fastq_files/${line}_1.fastq ~/Desktop/RNAseq_analysis/sample_files/fastq_files/${line}_2.fastq \
       ${paired}/${line}_R1_P.fastq ${unpaired}/${line}_R1_UP.fastq \
       ${paired}/${line}_R2_P.fastq ${unpaired}/${line}_R2_UP.fastq \
       ILLUMINACLIP:TruSeq4-PE.fa:2:10:5 \
      LEADING:30 TRAILING:30 SLIDINGWINDOW:2:30 MINLEN:60"
      
    trimmomatic PE -threads 6  ~/Desktop/RNAseq_analysis/sample_files/fastq_files/${line}_1.fastq \
                  ~/Desktop/RNAseq_analysis/sample_files/fastq_files/${line}_2.fastq \
       	  ${paired}/${line}_R1_P.fastq ${unpaired}/${line}_R1_UP.fastq \
                 ${paired}/${line}_R2_P.fastq ${unpaired}/${line}_R2_UP.fastq \
                 ILLUMINACLIP:TruSeq4-PE.fa:2:30:10:1:TRUE \
                 LEADING:30 TRAILING:30 SLIDINGWINDOW:2:30 MINLEN:45
                
                 
                    
     
     echo -e "\n"
       
 echo "**************************************  ${line} Trimmometry process completed **********************************************"
 echo -e "\n"
 
    
    
    
   
 
 
 
 
  echo "*******************************************Aligning ${line} with genome index*************************************************************************"
  echo -e "\n"
  echo "${n}"
  echo "bwa mem -t 5 ${index_path} ${paired}/${line}_R1_P.fastq ${paired}/${line}_R2_P.fastq -o ${aligned_file}/${line}.sam"
  bwa mem -t 5 ${index_path} ${paired}/${line}_R1_P.fastq ${paired}/${line}_R2_P.fastq -o ${aligned_file}/${line}.sam
  
  echo -e "\n"
  echo "*****************************Aligned finished for the file ${line}**********************************************************"
  echo -e "\n"
  
  aligned_file="/home/user/Desktop/RNAseq_analysis/sample_files/alignmed_Output"
  
  
  
  #features count 
  
  featureCounts -F "GTF" -t "gene" -a ${annotate_file} -o ${count_output}/${line}.txt  ${aligned_file}/${line}.sam 
 
 

 
steps_counts=$((steps_counts +1))
  
    

done<"$file"

echo " Total files ${steps_counts}"



