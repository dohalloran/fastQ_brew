#
# module for fastQ_brew
#
# Copyright Damien O'Halloran
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

fastQ_brew - a modern module for preprocessing of fastQ formatted files

=head1 DESCRIPTION

Provides methods for filtering and trimming reads by length and quality.

=head1 FEEDBACK

damienoh@gwu.edu

=head2 Mailing Lists

User feedback is an integral part of the evolution of this module. Send your comments and suggestions preferably to one of the mailing lists. Your participation is much appreciated.

=head2 Support

Please direct usage questions or support issues to:
<damienoh@gwu.edu>
Please include a thorough description of the problem with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the GitHub bug tracking system to help keep track of the bugs and their resolution.  Bug reports can be submitted via the GitHub page:

 https://github.com/dohalloran/fastQ_brew/issues

=head1 AUTHORS - Damien OHalloran

Email: damienoh@gwu.edu

=head1 APPENDIX

The rest of the documentation details each method

=cut

# Let the code begin...

package fastQ_brew;

use Moose;
use Modern::Perl;
with 'MooseX::Getopt', 'MooseX::Getopt::Usage::Role::Man';
use fastQ_brew_Utilities;
use fastQ_brew_Conversions;
use File::Copy qw(copy);
use Text::LevenshteinXS qw(distance);
use if ( $^O eq 'MSWin32' ), 'Win32::Console::ANSI';
use Term::ANSIColor qw(:constants);
use List::MoreUtils qw(first_index);
use List::MoreUtils qw{any};
use Sys::Hostname;
use POSIX qw/floor/;
use autodie;

##################################
our $VERSION = '2.0';
##################################

=head2 fastQ_brew->new_with_options()

 Function: Populates the user data into $self 
 Returns : nothing returned
 Args    :
  --i, input file (required)
  --x, input file2 only to check demultiplexing 
  --plex, checks that two FASTQ files were demultiplexed correctly
  --o, output file (required)
  --lib, library type  (default is sanger)       
  --dup, remove duplicate reads
  --qf, filter by phred (suggested default=20, min Q score=8 by default)
  --prob, probability that the read contains 0 errors (suggested default=0.5)
  --lf, filter by read length (suggested default =25)
  --trim_l, trim reads starting at left end
  --trim_r, trim reads starting at left end
  --truseq, removes truseq adapters from reads
  --adpt, remove a user supplied adapter     
  --fasta, convert to fastA format 
  --dna_rna, convert reads to RNA
  --rev_comp, reverse complement reads
  --no_n, remove non-designated bases from reads
  --help, Print this help

=cut

##################################

has 'i'        	=> ( is => 'rw', isa => 'Str',  required => 1 );
has 'x'         => ( is => 'rw', isa => 'Str',  required => 0 );
has 'plex'      => ( is => 'rw', isa => 'Bool', default => 0 );
has 'o'        	=> ( is => 'rw', isa => 'Str',  required => 1 );
has 'lib'      	=> ( is => 'rw', isa => 'Str',  default  => "sanger" );
has 'dup'      	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'qf'       	=> ( is => 'rw', isa => 'Num',  default  => 0 );
has 'prob'      => ( is => 'rw', isa => 'Num',  default  => 1 );
has 'lf'       	=> ( is => 'rw', isa => 'Int',  default  => 0 );
has 'truseq'   	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'trim_l'   	=> ( is => 'rw', isa => 'Int',  default  => 0 );
has 'trim_r'   	=> ( is => 'rw', isa => 'Int',  default  => 0 );
has 'adpt'     	=> ( is => 'rw', isa => 'Str',  default  => 0 );
has 'fasta'    	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'dna_rna'  	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'rev_comp' 	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'no_n'     	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'start'    	=> ( is => 'ro', isa => 'Int',  default  => time );
has help       	=> (
    is            	=> 'ro',
    isa           	=> 'Bool',
    default       	=> 0,
    documentation => qq{view driver script for input flag details}
);

##################################

=head2 truseq adapters

>TruSeq_Universal_Adapter
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq_Adapter_Index_1_6
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_2
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_3
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_4
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_5
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_6
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_7
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_8
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_9
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_10
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_11
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_12
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_13
GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_14
GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_15
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_16
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_18_7
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_19
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_20
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_21
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_22
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_23
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_25
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_27
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG

=cut

##################################

my %truseq_seqs = (
    0  => "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    1  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
    2  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
    3  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
    4  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
    5  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
    6  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
    7  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
    8  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
    9  => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
    10 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
    11 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
    12 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
    13 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG",
    14 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG",
    15 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG",
    16 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG",
    17 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG",
    18 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG",
    19 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG",
    20 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG",
    21 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG",
    22 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG",
    23 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG",
    24 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG"
);

##################################

=head2 run_fastQ_brew()

 Title   : run_fastQ_brew()
 Usage   : $self->run_fastQ_brew(%arg)
 Function: option to filter reads
 Returns : filtered fastQ file
 Args    : $self, %arg

=cut

##################################

sub run_fastQ_brew {
    my ( $self, %arg ) = @_;

    my ( $temp, @temp, $right_seq, $right_qual, $right_seqA, $right_qualA,
        @number_matchesR );
    my ( $left_seq, $left_qual, $left_seqA, $left_qualA, $count,
        @number_matchesL );

    # add user defined adapters
    unless ( $self->{adpt} eq 0 ) {
        $self->{adpt} =~ s/[^ACGTagct]//g;
        my $new_value = uc $self->{adpt};
        $truseq_seqs{25} = $new_value;
    }

    # log data on sys and infile
    my $inSize = -s $self->{i};
    print BOLD GREEN, "\nPerl Version: \t\t" . $^V;
    print "\nOperating system: \t" . $^O;
    print "\nHostname: \t\t" . hostname;
    print "\nfastQ_brew version: \t" . $VERSION;
    print "\nInput file size: \t" . $inSize . "KB", RESET;
    print BOLD CYAN, "\n\nFiltering file with fastQ_brew...", RESET;

    #check for file conversion calls or demultiplex check
    unless ( $self->{fasta} eq 0 ) {
        _convert_fasta( $self->{i} );
    }
    unless ( $self->{rev_comp} eq 0 ) {
        _rev_comp( $self->{i} );
    }
    unless ( $self->{dna_rna} eq 0 ) {
        _dna_rna( $self->{i} );
    }
    unless ( $self->{plex} eq 0 ) {
        _de_plex( $self->{i}, $self->{x} );
    }

    my $trimmed = "temp_";

    # new outfile resulting from filtering
    open my $fh_out, '>', $trimmed
      or die "Cannot open $trimmed: $!";

    open my $fh, '<', $self->{i}
      or die "Cannot open $self->{in_file}: $!";

    my %seen;
    my $unique;
    while (<$fh>) {
        $count++;
        chomp( $temp[0] = $_ );
        chomp( $temp[1] = <$fh> );
        chomp( $temp[2] = <$fh> );
        chomp( $temp[3] = <$fh> );

        # see if read is unique
        if ( $self->{dup} ne 0 && $seen{ $temp[1] }++ ) {
            $unique = 0;
        }
        else {
            $unique = 1;
        }

        my $key;

        ##  start by considering right end adapters
        if ( $self->{truseq} ne 0 || $self->{adpt} ne 0 ) {
            @number_matchesR = ();
            foreach $key ( sort { $a <=> $b } keys %truseq_seqs ) {
                my $hammingR = substr $temp[1], -length $truseq_seqs{$key},
                  length $truseq_seqs{$key};
                push @number_matchesR,
                  distance( $hammingR, $truseq_seqs{$key} );
            }

            if (   any { $_ == 0 } @number_matchesR
                or any { $_ == 1 } @number_matchesR )
            {
                my $index = first_index { $_ < 2 } @number_matchesR;
                $right_seqA  = substr $temp[1], 0, -length $truseq_seqs{$index};
                $right_qualA = substr $temp[3], 0, -length $truseq_seqs{$index};
            }
            else {
                $right_seqA  = $temp[1];
                $right_qualA = $temp[3];
            }
        }
        else {
            $right_seqA  = $temp[1];
            $right_qualA = $temp[3];
        }

        ##  left end adapters
        if ( $self->{truseq} ne 0 || $self->{adpt} ne 0 ) {
            @number_matchesL = ();
            foreach $key ( sort { $a <=> $b } keys %truseq_seqs ) {
                my $hammingL = substr $temp[1], 0, length $truseq_seqs{$key};
                push @number_matchesL,
                  distance( $hammingL, $truseq_seqs{$key} );
            }

            if (   any { $_ == 0 } @number_matchesL
                or any { $_ == 1 } @number_matchesL )
            {
                my $index = first_index { $_ < 2 } @number_matchesL;
                $left_seqA = substr $right_seqA, length $truseq_seqs{$index},
                  length $right_seqA;
                $left_qualA = substr $right_qualA, length $truseq_seqs{$index},
                  length $right_qualA;
            }
            else {
                $left_seqA  = $right_seqA;
                $left_qualA = $right_qualA;
            }
        }
        else {
            $left_seqA  = $right_seqA;
            $left_qualA = $right_qualA;
        }

        # if trim_l option not equal to N, then trim
        if ( $self->{trim_l} ne 0 ) {
            $left_seq = substr $right_seqA, $self->{trim_l}, length $right_seqA;
            $left_qual = substr $right_qualA, $self->{trim_l},
              length $right_qualA;

        }
        elsif ( $self->{trim_l} eq 0 ) {
            $left_seq  = $left_seqA;
            $left_qual = $left_qualA;
        }

        # if trim_r option not equal to N, then trim
        if ( $self->{trim_r} ne 0 ) {
            $right_seq  = substr $left_seq,  0, -$self->{trim_r};
            $right_qual = substr $left_qual, 0, -$self->{trim_r};

        }
        elsif ( $self->{trim_r} eq 0 ) {
            $right_seq  = $left_seq;
            $right_qual = $left_qual;
        }

        my $temp_prob;

        # get the read phred score
        if ( $self->{qf} ne 0 || $self->{prob} ne 0 ) {
            $temp_prob = _prob_calc( $right_qual, $self->{lib}, $self->{qf}, $self->{prob} );
        }
        else {
            $temp_prob = 1;
        }

        # conditionals for length, quality, de-dupe, and removing reads with N's
        if (
            (
                   ( length $right_seq > $self->{lf} && $self->{lf} ne 0 )
                || ( $self->{lf} eq 0 )
            )
            && $temp_prob == 1
            && (   ( $unique eq 1 && $self->{dup} ne 0 )
                || ( $self->{dup} eq 0 ) )
            && (   ( $right_seq !~ m/N/i && $self->{no_n} ne 0 )
                || ( $self->{no_n} eq 0 ) )
          )
        {
            print $fh_out "$temp[0]\n";
            print $fh_out "$right_seq\n";
            print $fh_out "$temp[2]\n";
            print $fh_out "$right_qual\n";
        }

    }
    $self->{original_reads} = $count;
    $self->{i}              = $trimmed;
    close $fh;
    close $fh_out;
    $self->_cleanup(%arg);
}

##################################

=head2 _cleanup()

 Title   : _cleanup()
 Usage   : _cleanup(%arg);
 Function: delete tmp files
 Returns : nothing
 Args    : %arg

=cut

##################################

sub _cleanup {
    my ( $self, %arg ) = @_;
    my $count;
    my $clean_out = $self->{o};
    copy $self->{i}, $clean_out;
    my $size = -s $clean_out;

    open my $fh, '<', $clean_out
      or die "Cannot open $clean_out: $!";
    print BOLD YELLOW,
      "\n\nInput number of reads: \t" . $self->{original_reads};
    print "\nOutput file size: \t" . $size . "KB";
    unless ( $size == 0 ) {
        $count += tr/\n/\n/ while sysread( $fh, $_, 2**16 );
        my $outCount = $count / 4;
        print "\nOutput number of reads:\t" . $outCount;
        print "\nReads removed:\t\t"
          . ( ( $self->{original_reads} ) - ($outCount) );
        print "\nPercent removed:\t"
          . floor(((( ( $self->{original_reads} ) - ($outCount) ))/$self->{original_reads})*100)."%";
    }
    my $duration = time - $self->{start};
    print "\nExecution time:\t\t" . $duration . "secs\n\n", RESET;

    # delete tmp files
    if ( -e "temp_" ) {
        unlink "temp_";
    }
    print BOLD RED, "fastQ_brew is now finished\n\n", RESET;
}

###################################

=head1 LICENSE AND COPYRIGHT

 Copyright (C) 2017 Damien M. O'Halloran
 GNU GENERAL PUBLIC LICENSE

=cut

1;
