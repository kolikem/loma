# loma

LoMA is a localized assembly tool for long reads.  
The current version is optimized for ONT reads although it can also be used for PacBio data.  
The following instructions assume users use UNIX-like OS. LoMA needs minimap2 (Heng Li) and MAFFT (Katoh et al.), so users please install them beforehand.

## Users' Guide
You can download the source code.
```sh
$ git clone https://github.com/kolikem/loma
```
Usage
LoMA runs with the following parameters:  
-I: input directory. LoMA takes FASTQ files as input data. Please make an input directory inclusing FASTQ files for which you want to obtain consensus sequences. Only one file in the directory works.  
-O: output directory. LoMA generates a consensus sequence for each region in input FASTQ files. The final sequences will be stored in a directory named "CONSENSUS", which is automatically made in the process.  
-b: block size. This is the window size. When LoMA makes a consensus sequene, it divides the whole region into pieces. (default=3000)  
-s: step size. This is the step size of blocks. (default=2000)  
-h: The number of reads for both ends' truncation. Both ends of a region typically have a lower depth, so LoMA limits the area constructed. It can be changed according to the size of users' data. (default=10)  
-d: The number of sigma in read classification. Users usually do not need to change this parameter. (default=3)


```sh
$ cd loma/loma
$ sh loma.sh <INPUT DIR> <OUTPUT DIR> 3 0 ont <minimap2> <MAFFT>
```
INPUT DIR is a directory in which FASTQ files that you want to assembly are stored.
Input format is FASTQ.
OUTPUT DIR is a directory in which you want to save consensus sequences.
For the last two parameters, please paste pathes to minimap2 and MAFFT on your computer.

## Dependency
minimap2 ver.2.0 <=  
MAFFT ver.7 <=
