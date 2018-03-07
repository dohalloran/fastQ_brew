#!/usr/bin/perl

use Moose;
use Modern::Perl;
use autodie;
use fastQ_brew;

my $app = fastQ_brew->new_with_options();
$app->run_fastQ_brew();

=head1 OPTIONS

  --i, input file (required)
  --x, input file2, only to check demultiplexing 
  --plex, checks that two FASTQ files were demultiplexed correctly
  --o, output file (required)    
  --lib, library type  (default is sanger)    
  --qf, filter by phred (suggested default=20, min Q score=8 by default and prob = 0.5)
  --lf, filter by read length (suggested default=25)
  --truseq, removes truseq adapters from reads
  --trim_l, trim reads starting at left end
  --trim_r, trim reads starting at left end   
  --fasta, convert to fastA format 
  --dna_rna, convert reads to RNA
  --rev_comp, reverse complement reads
  --help, Print this help
  
=cut
