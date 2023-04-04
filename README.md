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
fasta file
### Usage
A user can download all sources by git clone.
```sh
$ git clone https://github.com/kolikem/loma
```
Then you execute SETUP as follows.
```sh
$ cd loma
$ sh SETUP.sh
```
Now you are ready to use the tool. For general usage for reconstructions of localized genomic regions, a user can run LoMA by:
```sh
$ loma -I <INPUT> -O <OUTPUT>
```
INPUT is a directory designated by a user, which is supposed to include fastq file(s) from localized regions.
OUTPUT is also a directory defined by a user and will have three directories newly made: CONSENSUS, dir1, dir2. Final CSs are put in CONSENSUS directory with extension of .cs (fasta).
### Parameters
LoMA runs with the following parameters:  
  
-I: <PATH> input directory. LoMA takes FASTQ files as input data. Please make an input directory inclusing FASTQ files for which you want to obtain consensus sequences. Only one file in the directory works.  
  
-O: <PATH> output directory. LoMA generates a consensus sequence for each region in input FASTQ files. The final sequences will be stored in a directory named "CONSENSUS", which is automatically made in the process.  
  
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
minimap2 ver.2.0 <=  
MAFFT ver.7 <=
numpy

## Reference
Ikemoto, K., Fujimoto, H. & Fujimoto, A. Localized assembly for long reads enables genome-wide analysis of repetitive regions at single-base resolution in human genomes. Hum Genomics 17, 21 (2023). https://doi.org/10.1186/s40246-023-00467-7

