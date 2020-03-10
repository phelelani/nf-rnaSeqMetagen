#!/bin/bash

ftp_files=(ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/001/SRR5074531/SRR5074531_1.fastq.gz
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/001/SRR5074531/SRR5074531_2.fastq.gz 
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/000/SRR5074530/SRR5074530_1.fastq.gz
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/000/SRR5074530/SRR5074530_2.fastq.gz
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/009/SRR5074529/SRR5074529_1.fastq.gz
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/009/SRR5074529/SRR5074529_2.fastq.gz
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/008/SRR5074528/SRR5074528_1.fastq.gz
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/008/SRR5074528/SRR5074528_2.fastq.gz)

for url in "${ftp_files[@]}"
do
    lftp -e 'pget -n 5 -c '"${url}"'; bye'
done
