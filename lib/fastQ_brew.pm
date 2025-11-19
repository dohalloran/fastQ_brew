package fastQ_brew;

use strict;
use warnings;
use feature qw(say);

use Moose;

use fastQ_brew_Utilities;
use fastQ_brew_Conversions;
use File::Copy qw(copy);
use Sys::Hostname;
use Term::ANSIColor qw(:constants);
use autodie;

##################################
our $VERSION = '2.1_core_only';
##################################

=head1 NAME

fastQ_brew - preprocessing of fastQ formatted files (core-only, parallel-safe)

=head1 DESCRIPTION

Filters and trims reads by length and quality, with optional TruSeq adapter
trimming and simple format conversions. This version is designed to run with
only core Perl plus Moose and Term::ANSIColor, and to be safe when used from
an external parallel wrapper.

=cut

#----------------------#
#   Attributes / CLI   #
#----------------------#

has 'i'        	=> ( is => 'rw', isa => 'Str',  required => 1 );  # input FASTQ
has 'x'         => ( is => 'rw', isa => 'Str',  required => 0 );  # for plex
has 'plex'      => ( is => 'rw', isa => 'Bool', default  => 0 );
has 'o'        	=> ( is => 'rw', isa => 'Str',  default  => "filtered.fq" );
has 'lib'      	=> ( is => 'rw', isa => 'Str',  default  => "sanger" );
has 'qf'       	=> ( is => 'rw', isa => 'Num',  default  => 0 );  # error prob cutoff (%)
has 'lf'       	=> ( is => 'rw', isa => 'Int',  default  => 0 );  # length filter
has 'truseq'   	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'trim_l'   	=> ( is => 'rw', isa => 'Int',  default  => 0 );
has 'trim_r'   	=> ( is => 'rw', isa => 'Int',  default  => 0 );
has 'fasta'    	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'dna_rna'  	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'rev_comp' 	=> ( is => 'rw', isa => 'Bool', default  => 0 );
has 'start'    	=> ( is => 'ro', isa => 'Int',  default  => time );

# keep workers for backwards compatibility (not used in this core-only build)
has 'workers'   => (
    is            => 'rw',
    isa           => 'Int',
    default       => 1,
    documentation => 'ignored in core-only build; external parallel wrapper is recommended',
);

# used internally to track the temporary output file
has 'tmp_file'  => (
    is        => 'rw',
    isa       => 'Str',
    predicate => 'has_tmp_file',
);

has 'help'     	=> (
    is            	=> 'ro',
    isa           	=> 'Bool',
    default       	=> 0,
    documentation => 'see fastQ_brew.pl for CLI usage',
);

#----------------------#
#    TruSeq adapters   #
#----------------------#

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
    24 => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG",
);

#----------------------#
#   Helper functions   #
#----------------------#

# simple Hamming-like distance (length difference counts as mismatches)
sub _simple_distance {
    my ( $a, $b ) = @_;
    my $la = length $a;
    my $lb = length $b;
    my $len = $la < $lb ? $la : $lb;
    my $dist = 0;

    for ( my $i = 0; $i < $len; $i++ ) {
        my $ca = substr( $a, $i, 1 );
        my $cb = substr( $b, $i, 1 );
        $dist++ if $ca ne $cb;
    }
    $dist += abs( $la - $lb );
    return $dist;
}

# process one FASTQ record: [header, seq, plus, qual]
sub _process_record {
    my ( $self, $rec_ref ) = @_;
    my @temp = @$rec_ref;

    # TruSeq adapter removal
    if ( $self->{truseq} ) {
        my @number_matchesR;
        my @number_matchesL;

        foreach my $key ( sort { $a <=> $b } keys %truseq_seqs ) {
            my $adapter = $truseq_seqs{$key};

            my $hammingR = substr $temp[1], -length($adapter), length($adapter);
            push @number_matchesR, _simple_distance( $hammingR, $adapter );

            my $hammingL = substr $temp[1], 0, length($adapter);
            push @number_matchesL, _simple_distance( $hammingL, $adapter );
        }

        # right side
        my $indexR = -1;
        for ( my $i = 0 ; $i < @number_matchesR ; $i++ ) {
            if ( $number_matchesR[$i] < 2 ) {
                $indexR = $i;
                last;
            }
        }
        if ( $indexR != -1 ) {
            my $adapter = $truseq_seqs{$indexR};
            $temp[1]  = substr $temp[1], 0, -length($adapter);
            $temp[3] = substr $temp[3], 0, -length($adapter);
        }

        # left side
        my $indexL = -1;
        for ( my $i = 0 ; $i < @number_matchesL ; $i++ ) {
            if ( $number_matchesL[$i] < 2 ) {
                $indexL = $i;
                last;
            }
        }
        if ( $indexL != -1 ) {
            my $adapter = $truseq_seqs{$indexL};
            $temp[1] = substr $temp[1], length($adapter), length $temp[1];
            $temp[3] = substr $temp[3], length($adapter), length $temp[3];
        }
    }

    # Left trim
    if ( $self->{trim_l} ) {
        my $tl = $self->{trim_l};
        $temp[1] = substr $temp[1], $tl, length $temp[1];
        $temp[3] = substr $temp[3], $tl, length $temp[3];
    }

    # Right trim
    if ( $self->{trim_r} ) {
        my $tr = $self->{trim_r};
        $temp[1]  = substr $temp[1], 0, -$tr;
        $temp[3] = substr $temp[3], 0, -$tr;
    }

    # Quality filter
    my $temp_prob = 1;
    if ( $self->{qf} ) {
        $temp_prob = _prob_calc( $temp[3], $self->{lib}, $self->{qf} );
    }

    # Length + prob condition
    if (
        (
               ( $self->{lf} && length( $temp[1] ) > $self->{lf} )
            || ( !$self->{lf} )
        )
        && $temp_prob == 1
      )
    {
        return \@temp;
    }

    return;    # filtered out
}

#----------------------#
#        API          #
#----------------------#

sub run_fastQ_brew {
    my ( $self, %arg ) = @_;

    my $inSize = -s $self->{i};
    print BOLD GREEN, "\nPerl Version: \t\t" . $^V;
    print "\nOperating system: \t" . $^O;
    print "\nHostname: \t\t" . hostname;
    print "\nfastQ_brew version: \t" . $VERSION;
    print "\nInput file size: \t" . $inSize . " bytes", RESET;
    print BOLD CYAN, "\n\nFiltering file with fastQ_brew...", RESET;

    # conversions / demultiplex check
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

    # core-only build: serial path only (external parallel is recommended)
    $self->_run_serial(%arg);
}

sub _run_serial {
    my ( $self, %arg ) = @_;

    # UNIQUE temp file per process (PID-based)
    my $trimmed = "temp_$$";

    my @temp;
    my $count = 0;

    open my $fh_out, '>', $trimmed
      or die "Cannot open $trimmed: $!";

    open my $fh, '<', $self->{i}
      or die "Cannot open $self->{i}: $!";

    while (1) {
        my $h = <$fh>;
        my $s = <$fh>;
        my $p = <$fh>;
        my $q = <$fh>;
        last unless defined $q;   # guard against truncated FASTQ

        chomp( $temp[0] = $h );
        chomp( $temp[1] = $s );
        chomp( $temp[2] = $p );
        chomp( $temp[3] = $q );

        $count++;

        my $out_rec = $self->_process_record( \@temp );
        next unless $out_rec;

        print $fh_out "$out_rec->[0]\n";
        print $fh_out "$out_rec->[1]\n";
        print $fh_out "$out_rec->[2]\n";
        print $fh_out "$out_rec->[3]\n";
    }

    $self->{original_reads} = $count;
    $self->{i}              = $trimmed;
    $self->{tmp_file}       = $trimmed;

    close $fh;
    close $fh_out;

    $self->_cleanup(%arg);
}

sub _cleanup {
    my ( $self, %arg ) = @_;
    my $count = 0;

    my $clean_out = $self->{o};
    copy $self->{i}, $clean_out;

    my $size = -s $clean_out;

    open my $fh, '<', $clean_out
      or die "Cannot open $clean_out: $!";

    print BOLD YELLOW,
      "\n\nInput number of reads: \t" . $self->{original_reads};
    print "\nOutput file size: \t" . $size . " bytes";

    unless ( $size == 0 ) {
        $count += tr/\n/\n/ while sysread( $fh, $_, 2**16 );
        my $outCount = $count / 4;
        print "\nOutput number of reads:\t" . $outCount;
        print "\nReads removed:\t\t"
          . ( ( $self->{original_reads} ) - ($outCount) );
        print "\nPercent removed:\t"
          . ( sprintf "%.6f",
            ( ( ( $self->{original_reads} - $outCount ) / $self->{original_reads} ) * 100 ) ) . "%";
    }
    my $duration = time - $self->{start};
    print "\nExecution time:\t\t" . $duration . " secs\n\n", RESET;

    # delete this instance's temp file
    if ( $self->has_tmp_file && -e $self->{tmp_file} ) {
        unlink $self->{tmp_file};
    }

    print BOLD RED, "fastQ_brew is now finished\n\n", RESET;
}

1;
