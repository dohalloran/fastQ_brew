# `fastQ_brew v.2.0`

[![GitHub license](https://img.shields.io/badge/license-GPL_2.0-orange.svg)](https://raw.githubusercontent.com/dohalloran/fastQ_brew/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/dohalloran/fastQ_brew.svg)](https://github.com/dohalloran/fastQ_brew/issues)

- [x] `Pre-processing of FASTQ [1] reads`
- [x] `Check that files were demultiplxed correctly`
- [x] `Filter reads by length `
- [x] `Filter reads by quality`
- [x] `Trim reads`
- [x] `Removes standard Truseq adapters`
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
  #standard use
  perl brew_driver.pl -i <input_file> -o <output_file> --qf 20 --lf 25 --truseq
  
  #check that files were demultiplxed correctly
  perl brew_driverl.pl -i <input_file1> -x <input_file2> -o <output_file> --plex
  
  #see below for command flags 
```

## Command Line Arguments
### Filtering Options
 ```perl   
#A FASTQ read is removed if the following criteria are not met for a given read: 
# i) average read Phred Quality [2-3] is not above the user supplied threshold (see --qf, suggested default=20)
# ii)the min Phred Q score for any given nucleotide is not above 8
# iii)the probability that the read contains 0 errors i.e. E < 1 (default=0.5)
#only reads with an average Q score above --qf will be kept
        --qf 20
#filter by read length - reads below this length will be removed       
        --lf 25
#remove x bases from left end of every read 
        --trim_l 5
#remove x bases from right end of every read
        --trim_r 3
#remove standard truseq adapters (permits 1 mismatch) from both ends (very slow!)
        --truseq
```

### File Conversions and de-multiplex check
 ```perl   
#check that 2 FASTQ files were demultiplexed correctly 
#fastQ_brew outputs the barcodes for each file and compares (union and intersection) between two files 
        --plex
        -i <input_file1>
        -x <input_file2>
        -o <output_file>
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
        --i reads.fastq
#output FASTQ file (required) 
        --o filtered.fq
#library type i.e. sanger (default) or illumina 
        --lib sanger
#print flag options to stdout
        --help  
```

### References
1. Cock PJ, Fields CJ, Goto N, Heuer ML, Rice PM. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Res. 2010;38(6):1767–71

2. Ewing B, Hillier L, Wendl MC, Green P. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Res. 1998;8(3):175–85.

3. Ewing B, Green P. Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Res. 1998;8(3):186–94.

## Contributing
All contributions are welcome.

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/fastQ_brew/issues).

## License 
GNU GENERAL PUBLIC LICENSE





