#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use lib 'lib';
use fastQ_brew;

# default options (matching Moose defaults where sensible)
my %opt = (
    o       => 'filtered.fq',
    lib     => 'sanger',
    qf      => 0,
    lf      => 0,
    truseq  => 0,
    trim_l  => 0,
    trim_r  => 0,
    fasta   => 0,
    dna_rna => 0,
    rev_comp=> 0,
    plex    => 0,
);

GetOptions(
    'i=s'       => \$opt{i},       # input FASTQ (required)
    'x=s'       => \$opt{x},       # second input (for plex)
    'o=s'       => \$opt{o},       # output file
    'lib=s'     => \$opt{lib},     # sanger / illumina
    'qf=f'      => \$opt{qf},      # quality cutoff (percent error probability)
    'lf=i'      => \$opt{lf},      # length cutoff
    'truseq!'   => \$opt{truseq},  # remove TruSeq adapters
    'trim_l=i'  => \$opt{trim_l},  # left trim
    'trim_r=i'  => \$opt{trim_r},  # right trim
    'fasta!'    => \$opt{fasta},   # convert to FASTA
    'dna_rna!'  => \$opt{dna_rna}, # DNA->RNA
    'rev_comp!' => \$opt{rev_comp},# reverse complement
    'plex!'     => \$opt{plex},    # demultiplex check
) or die "Error in command line arguments\n";

die "Input FASTQ required (--i)\n"
    unless defined $opt{i};

# Don't pass x at all unless it was provided
delete $opt{x} unless defined $opt{x};

my $app = fastQ_brew->new(%opt);
$app->run_fastQ_brew();
