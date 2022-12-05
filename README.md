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
You can use LoMA by running the shell script.
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
