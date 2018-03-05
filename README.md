# `fastQ_brew v.2.0`

[![GitHub license](https://img.shields.io/badge/license-GPL_2.0-orange.svg)](https://raw.githubusercontent.com/dohalloran/fastQ_brew/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/dohalloran/fastQ_brew.svg)](https://github.com/dohalloran/fastQ_brew/issues)

- [x] `Pre-processing of FASTQ [1] reads` 
- [x] `Filters reads by length `
- [x] `Filters reads by quality`
- [x] `Trims read ends`
- [x] `Removes duplicate reads`
- [x] `Removes standard Truseq adapters`
- [x] `Removes user supplied adapter`
- [x] `Performs various file conversions` 

![fastQ_brew LOGO](https://cloud.githubusercontent.com/assets/8477977/22077145/f29a177e-dd80-11e6-86a6-a211e8e1e103.jpg)

## Installation
1. Download and extract the fastQ_brew.zip file  
`tar -xzvf fastQ_brew.zip` 
or 
`git clone https://github.com/dohalloran/fastQ_brew.git`
2. The extracted dir will be called fastQ_brew  
  `cd fastQ_brew`   
  `perl Makefile.PL`  
  `make`  
  `make test`  
  `make install`  

## Usage 
### Type the following to run:  
 ```bash 
  #brew_driver.pl is a driver script within the lib folder 
  
  perl brew_driver.pl --i=<input_file> --o=<output_file> --qf=10 --lf=25 --truseq 
  
  #see below for command flags 
```

## Command Line Arguments
### Filtering Options
 ```perl   
#A FASTQ read is removed if the following two criteria are not both met for a given read: 
# i) average read Phred Quality [2-3] is not above the user supplied threshold (see --qf flag)
# ii)the min Phred is not within a specified number of standard deviations of the mean (see --dist flag)
#Quality filtering is based on the 68-95-99.7 rule [4] (equation below) with 1 std deviation should account
#for ~68% of data if normally distributed (which it's not, see below)
        --qf=10
#number of standard deviations between min and mean Phred score permitted (default is 5; lower for stringency)
        --dist=5
#filter by read length - reads below this length will be removed       
        --lf=25
#remove x bases from left end of every read 
        --trim_l=5
#remove x bases from right end of every read
        --trim_r=3
#remove user specified adapter (permits 1 mismatch) from both ends
        --adpt=AATGATACGGCGACCACCGAGATCTACACT
#remove standard truseq adapters (permits 1 mismatch) from both ends
        --truseq
#remove duplicate reads 
        --dup
#remove reads that contain non designated bases e.g. N or [^ATGC] 
        --no_n
```
 <a href="https://www.codecogs.com/eqnedit.php?latex=Pr(\mu&space;-&space;\sigma&space;\leq&space;X&space;\leq&space;\mu&space;&plus;&space;\sigma&space;)&space;\approx&space;0.6827" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Pr(\mu&space;-&space;\sigma&space;\leq&space;X&space;\leq&space;\mu&space;&plus;&space;\sigma&space;)&space;\approx&space;0.6827" title="Pr(\mu - \sigma \leq X \leq \mu + \sigma ) \approx 0.6827"/></a>
 

### File Conversions
 ```perl   
#convert FASTQ file to FASTA format file
        --fasta
#convert the DNA for each read to RNA 
        --dna_rna
#reverse complement the FASTQ reads 
        --rev_comp
```

### Odds and Ends
 ```perl   
#input FASTQ file (required) 
        --i=reads.fastq
#output FASTQ file (required) 
        --o=filtered.fq
#library type i.e. sanger (default) or illumina 
        --lib=sanger
#print flag options to stdout
        --help  
```
### FASTQ Distribution

![fastQ brew histogram](https://user-images.githubusercontent.com/8477977/30244895-1e35c1bc-9596-11e7-9de8-b7c3c074f3e9.png)
![fastQ_brew dist](https://user-images.githubusercontent.com/8477977/30244876-bdcf8ce0-9595-11e7-9f0a-cecf82622ad3.png)

### References
1. Cock PJ, Fields CJ, Goto N, Heuer ML, Rice PM. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Res. 2010;38(6):1767–71

2. Ewing B, Hillier L, Wendl MC, Green P. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Res. 1998;8(3):175–85.

3. Ewing B, Green P. Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Res. 1998;8(3):186–94.

4. [Abraham de Moivre. The doctrine of chances: or, A method of calculating the probabilities of events in play, 1667-1754](https://archive.org/details/doctrineofchance00moiv)

## Contributing
All contributions are welcome.

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/fastQ_brew/issues).

## License 
GNU GENERAL PUBLIC LICENSE





