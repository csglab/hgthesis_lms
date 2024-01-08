#############
DESCRIPTION
#############

This directory is organized by dataset (sims = simulated data, tcga = TCGA data), and contains the scripts and outputs described in sections:
  - 2.1.1 Preparing TRex inputs
  - 2.1.2 TRex statistical model 

In addition, the `references` folder contains the all genome reference annotations required for this stage. 
At this stage, the inputs to the pipeline are raw fastq files. The raw fastqs are not part of the repository due to file size limitations, but the code used to process them with salmon is provided in the scripts folder. 

##############
FOLDER DETAILS
##############

Each dataset folder follows this basic structure:
    - output = folder containing all output files from multiple processing steps
    - scripts = code necessary to generate files in output

## sims ##

For coding purposes, the confounder variable in the simulated data is called batch. 
The outputs in this folder are organized by run. Each run corresponds to one of 5 random simulations (SIMULATION_NUMBER) performed for every value of confounder effect strength (STRENGTH). 

This information is captured by the name of the directory as:

run_[STRENGTH]_rs[SIMULATION_NUMBER]

STRENGTH = 0,0.25,0.5,0.75 or 1
SIMULATION_NUMBER = 1,2,3,4 or 5

Inside each run number there are three directories containing the output of each program: salmon, tximport and trex. The trex folder contains two RData files with the event counts in two different formats (tximport and tximeta). 

salmon
------
These results are structured by sample. Each sample contains the standard outputs of a salmon run.

tximport
--------
Contains one .RData file with the tximport object and one .csv files with the raw transcript counts.

trex
----

The results files from TRex (with the extension .tsv) are organized with the following naming convention: 

[EVENT].[FILE_TYPE].[MODEL_NAME].[CONTRAST_NAME].tsv

EVENT = SE
FILE_TYPE = res
MODEL_NAME = cell_line (no confounder added), cell_line.batch (with confounder)
CONTRAST_NAME = trex_countA, cell_lineSRSF9.KD, batch, trex_countA.cell_lineSRSF9.KD, or trex_countA.batch  

## tcga ## 

This folder is structured per cancer type. Each cancer type is a folder that contains two subdirectories with the outputs of each programs executed at this stage: tximport and trex.

Due to the large number of files, the raw salmon outputs are not provided here, but the output structure is the same as in the sims dataset. The outpus of each cancer type has the following structure:

tximport
--------

There are two .RData files containing tximport objects. The first one named [CANCER].tximport.RData corresponds to transcript-based counts and the second one named [CANCER].gene.tximport.RData contains the gene-level counts. Additionally, there is a .csv file containing the library sizes of all samples (named by their TCGA file id). 

trex
----

The trex output is organized into three subfolders: 
1. event_counts = contains a .csv and a .RData file with the event counts and the tximport-like object with event counts of each type of AS.
2. coefs_stage = has the resulting coefficients of the stage models. Individual file naming structure follows the same logic as the TRex results files in the sims folder.
3. coefs_condition = resulting coefficients of the condition models. Similar naming convention as 2. 

## scripts ##

The script files in each dataset folder operate as a pipeline controlled by the .sh files. The .sh files need to be executed according to the order indicated by the first two numbers of the file name. Internally, the pipeline files call the R scripts (.r files) and other unnumbered scripts located in the same directory. To understand what each of the unnumbered scripts do, look for the pipeline file in which it is being executed.


