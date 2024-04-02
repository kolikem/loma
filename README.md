# loma

## Summary
Overview:  
LoMA is a localized assembly tool for long reads.  
It starts with interpretations of the all-to-all alignment to a region of interest. LoMA lays out long reads after some filterations. It divides the layout into multiple blocks to make partial consensus sequences. They will be contatenated into one consensus sequence.

Note:  
The current version is optimized for ONT reads although it can also be used for PacBio data.  
The following instructions assume users of linux-like OS. Users who use Windows are recommended to use WSL to run commands.  
LoMA needs minimap2 (Heng Li) and MAFFT (Katoh et al.), so users please install them beforehand.  
Typically, LoMA takes 10-1000 reads.

## Guidance
### Input
fastq file
### Output
fasta file (.cs)
### Usage
Users can download all source files from "Releases". Please download the latest version v1.1.3.  
Then decompress and go to the directory loma.
```sh
$ cd loma
```
(case 1) You can run LoMA after executing SETUP.sh, which is for the establishment of the path.
```sh
$ sh SETUP.sh
```
Now you are ready to use the tool. For general usage for reconstructions of localized genomic regions, a user can run LoMA by:
```sh
$ loma -I <INPUT> -O <OUTPUT> -H <minimap2> -K <mafft>
```
(case 2) You can run LoMA without using SETUP.sh just by running the sh file.
```sh
$ sh loma -I <INPUT> -O <OUTPUT> -H <minimap2> -K <mafft>
```
INPUT is a directory designated by a user, which is supposed to include fastq file(s) from localized regions. Please make sure that INPUT is an absolute path. 
OUTPUT is also a directory defined by a user and will have three directories newly made; CONSENSUS, dir1, dir2. Final CSs are put in CONSENSUS directory with extension of .cs (fasta). Please make sure that OUTPUT is an absolute path as well.  
H and K are not necessary if their paths are reachable.
### Parameters
LoMA runs with the following parameters:  
  
-I: <PATH> input directory (absolute path). LoMA takes FASTQ files as input. Please make an input directory which contains at least one FASTQ file for which you want to obtain consensus sequence(s).   
  
-O: <PATH> output directory (absolute path). LoMA generates a consensus sequence for each input. The final sequences will be stored in a directory named "CONSENSUS", which will be made in the process.  
  
-b: <INT> block size. This is the window size. When LoMA makes a consensus sequene, it divides the whole region into pieces. (default=3000)  
  
-s: <INT> step size. This is the step size of blocks. (default=2000)  
  
-h: <INT> The number of reads for both ends' truncation. Both ends of a region typically have a lower depth, so LoMA limits the area constructed. It can be changed according to the size of users' data. (default=10)  
  
-d: <INT> The number of sigma in read classification. Users usually do not need to change this parameter. (default=3)  
  
-l: <ont/pb> Data. Nanopore (ONT) or PacBio. (default=ont)  
  
-c: <INT> Minimum coverage rate of a block. (default=0.7)  
  
-r: <FLOAT> A parameter for data filtering by alignment accuracy's rank (discard 100*x% of alignments). (default=0.5)  

-m: <INT> A parameter for data filtering by the number of matching bases. (default=1000)  
  
-H: <PATH> minimap2 (Heng Li). If you have not set the path.  
  
-K: <PATH> MAFFT (Katoh et al.). If you have not set the path.  

## Dependency
python >= 3.8  
minimap2 >= ver.2.0  
MAFFT >= ver.7  
numpy (python library)  
matplotlib (python library)  

## Reference
Ikemoto, K., Fujimoto, H. & Fujimoto, A. Localized assembly for long reads enables genome-wide analysis of repetitive regions at single-base resolution in human genomes. Hum Genomics 17, 21 (2023). https://doi.org/10.1186/s40246-023-00467-7

