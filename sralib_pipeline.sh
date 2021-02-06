#!/usr/bin/bash
# sra_pipeline.sh
# read list of SRA ids,
# get data using SRA toolkit
# process data and do simple stats using a python script


sra_list=$1


for id in `cat ${sra_list}`
do
    # Get SRA data
    # --------------------------------------------------
    # files are named: ${id}.fastq

    # call sra toolkit
    fastq-dump --split-spot --readids --dumpbase --skip-technical --clip --read-filter pass ${id}

    # Process fastq using python
    # --------------------------------------------------

    myprogramplot.py -i ${id}.fastq -o -g
    
done

# --------------------------------------------------
# alternative
# --------------------------------------------------

#while read -r id;
#do 
#    sra-toolkit blabla -input ${id};
#    MYPYTHINSCRIPT.py -par1 -par2 -par3 -input ${id}.fastq;
#done < ${sra_list}

