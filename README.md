# `fastQ_brew v.2.0`

[![GitHub license](https://img.shields.io/badge/license-GPL_2.0-orange.svg)](https://raw.githubusercontent.com/dohalloran/fastQ_brew/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/dohalloran/fastQ_brew.svg)](https://github.com/dohalloran/fastQ_brew/issues)

- [x] `Pre-processing of FASTQ reads`
- [x] `Check that files were demultiplexed correctly`
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
```cmd  
cd fastQ_brew   
perl Makefile.PL  
make  
make test  
make install  
```

## Usage 
### To run:  
```perl
use fastQ_brew;
use Moose;
use Modern::Perl;
use autodie;

my $app = fastQ_brew->new_with_options();
$app->run_fastQ_brew(); 
#see below for command flags
``` 
## Command Line Arguments
### Filtering Options
 ```perl   
#set the max probability that a fastQ [1] read will contain errors: suggested p<=50% (must be 1-100)
        --qf 50
#filter by read length - reads below this length will be removed       
        --lf 35
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
        --i <input_file>
#output FASTQ file (by default called filtered.fq) 
        --o <output_file>
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





