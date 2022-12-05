# loma

## Summary
Overview:  
LoMA is a localized assembly tool for long reads.  
It starts with interpretations of the all-to-all alignment to your ROI. LoMA lays out long reads after some filterations. It divides the layout into multiple blocks to make partial consensus sequences. They will be contatenated into one consensus sequence.

Note:  
The current version is optimized for ONT reads although it can also be used for PacBio data.  
The following instructions assume users use UNIX-like OS. LoMA needs minimap2 (Heng Li) and MAFFT (Katoh et al.), so users please install them beforehand.  
Typically, LoMA takes 10-1000 reads.

## Guidance
You can download the source code.
```sh
$ git clone https://github.com/kolikem/loma
```
For general usage, a user can run LoMA by:
```sh
$ sh loma.sh -I <INPUT> -O <OUTPUT>
```
LoMA runs with the following parameters:  
-I: <PATH> input directory. LoMA takes FASTQ files as input data. Please make an input directory inclusing FASTQ files for which you want to obtain consensus sequences. Only one file in the directory works.  
-O: <PATH> output directory. LoMA generates a consensus sequence for each region in input FASTQ files. The final sequences will be stored in a directory named "CONSENSUS", which is automatically made in the process.  
-b: <INT> block size. This is the window size. When LoMA makes a consensus sequene, it divides the whole region into pieces. (default=3000)  
-s: <INT> step size. This is the step size of blocks. (default=2000)  
-h: <INT> The number of reads for both ends' truncation. Both ends of a region typically have a lower depth, so LoMA limits the area constructed. It can be changed according to the size of users' data. (default=10)  
-d: <INT> The number of sigma in read classification. Users usually do not need to change this parameter. (default=3)  
-l: <ont/pb> Data. Nanopore (ONT) or PacBio. (default=ont)  
-c: <INT> Minimum coverage rate of a block. (default=0.7)  
-r: <FLOAT> A parameter for data filtering by alignment accuracy's rank %. (default=0.5)  
-m: <INT> A parameter for data filtering by the number of matching bases. (default=1000)  
-H: <PATH> minimap2 (Heng Li). If you have not set the path.  
-K: <PATH> MAFFT (Katoh et al.). If you have not set the path.  

INPUT DIR is a directory in which FASTQ files that you want to assembly are stored.
Input format is FASTQ.
OUTPUT DIR is a directory in which you want to save consensus sequences.
For the last two parameters, please paste pathes to minimap2 and MAFFT on your computer.

## Dependency
minimap2 ver.2.0 <=  
MAFFT ver.7 <=
