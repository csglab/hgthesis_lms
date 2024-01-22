#############
DESCRIPTION
#############

This directory contains all the necessary elements to conduct the benchmark of TRex against rMATS and SUPPA2 using custom synthetic datasets. It covers the sections: 
    - 2.2.1 Modeling confounding effects
    - 2.2.2 Computing performance metrics
    

##############
FOLDER DETAILS
##############

The basic structure is:
    - input = directory with input files required at the beginning of the pipeline (i.e samples metadata, fastq files from ENCODE (accessions: ENCSR972AZD, ENCSR341TTW,  and the rsem quantifications)
    - output = folder containing all output files from multiple processing steps. The name of the output folder describes the tool or tools used to generate the corresponding files. 
    - scripts = code necessary to generate files in output. The pipeline consists of serialized bash wrappers that run individual tools/scripts in the directory. The wrappers are named starting with the number of the corresponding step. To understand all the scripts in the folder follow the order of the wrappers. 