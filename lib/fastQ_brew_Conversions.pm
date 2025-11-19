package fastQ_brew_Conversions;

use strict;
use warnings;
use feature qw(say);

use Moose;
use base 'Exporter';
use Term::ANSIColor qw(:constants);
use autodie;

our @EXPORT = qw/ _convert_fasta _rev_comp _dna_rna _de_plex /;

#----------------------#
#  FASTQ -> FASTA      #
#----------------------#

sub _convert_fasta {
    my $in_file = shift or die "No input file provided to _convert_fasta\n";

    print BOLD CYAN, "\n\nConverting fastQ file to fastA...", RESET;

    my @temp;
    my $fasta = "fastA_convert.fa";

    open my $fh, '<', $in_file
      or die "Cannot open $in_file: $!";

    open my $fh_out, '>', $fasta
      or die "Cannot open $fasta: $!";

    while (<$fh>) {
        chomp( $temp[0] = $_ );
        chomp( $temp[1] = <$fh> );
        chomp( $temp[2] = <$fh> );
        chomp( $temp[3] = <$fh> );

        print $fh_out ">$temp[0]\n";
        print $fh_out "$temp[1]\n";
    }
    close $fh;
    close $fh_out;
}

#----------------------#
#   Reverse complement #
#----------------------#

sub _rev_comp {
    my $in_file = shift or die "No input file provided to _rev_comp\n";

    print BOLD CYAN, "\n\nReverse complementing fastQ reads...", RESET;

    my @temp;
    my $revComp;
    my $revcomp = "rev_comp.fastq";

    open my $fh, '<', $in_file
      or die "Cannot open $in_file: $!";

    open my $fh_out, '>', $revcomp
      or die "Cannot open $revcomp: $!";

    while (<$fh>) {
        chomp( $temp[0] = $_ );
        chomp( $temp[1] = <$fh> );
        chomp( $temp[2] = <$fh> );
        chomp( $temp[3] = <$fh> );

        $temp[1] =~ tr/ATGCatgc/TACGtacg/;
        $revComp = reverse( $temp[1] );

        print $fh_out "$temp[0]\n";
        print $fh_out "$revComp\n";
        print $fh_out "$temp[2]\n";
        print $fh_out "$temp[3]\n";
    }
    close $fh;
    close $fh_out;
}

#----------------------#
#   DNA -> RNA         #
#----------------------#

sub _dna_rna {
    my $in_file = shift or die "No input file provided to _dna_rna\n";

    print BOLD CYAN, "\n\nTranscribing fastQ reads...", RESET;

    my @temp;
    my $dna_to_rna = "dna_to_rna.fastq";

    open my $fh, '<', $in_file
      or die "Cannot open $in_file: $!";

    open my $fh_out, '>', $dna_to_rna
      or die "Cannot open $dna_to_rna: $!";

    while (<$fh>) {
        chomp( $temp[0] = $_ );
        chomp( $temp[1] = <$fh> );
        chomp( $temp[2] = <$fh> );
        chomp( $temp[3] = <$fh> );

        $temp[1] =~ tr/Tt/Uu/;

        print $fh_out "$temp[0]\n";
        print $fh_out "$temp[1]\n";
        print $fh_out "$temp[2]\n";
        print $fh_out "$temp[3]\n";
    }
    close $fh;
    close $fh_out;
}

#----------------------#
#   Demultiplex check  #
#----------------------#

sub _de_plex {
    my (@fileArray) = @_;
    my $in_file  = $fileArray[0] or die "No first input file provided to _de_plex\n";
    my $in_file2 = $fileArray[1] or die "No second input file provided to _de_plex\n";

    print BOLD CYAN, "\n\nChecking fastQ files for proper demultiplexing...", RESET;

    open my $fh, '<', $in_file
      or die "Cannot open $in_file: $!";

    open my $fh2, '<', $in_file2
      or die "Cannot open $in_file2: $!";

    my @tagsArr;
    my @tagsArr2;

    while ( my $line = <$fh> ) {
        if ( $line =~ m/:([NAGCT][NAGCT][NAGTC][NAGCT][NAGCT][NAGTC])$/ ) {
            push( @tagsArr, $1 ) unless grep { $_ eq $1 } @tagsArr;
        }
    }
    while ( my $line2 = <$fh2> ) {
        if ( $line2 =~ m/:([NAGCT][NAGCT][NAGTC][NAGCT][NAGCT][NAGTC])$/ ) {
            push( @tagsArr2, $1 ) unless grep { $_ eq $1 } @tagsArr2;
        }
    }

    if ( scalar @tagsArr == 0 ) {
        print BOLD CYAN, "\n\nNo tags found from file1\n", RESET;
    }
    else {
        print BOLD CYAN, "\n\nList of tags from file1: " . "@tagsArr\n", RESET;
    }

    if ( scalar @tagsArr2 == 0 ) {
        print BOLD CYAN, "\nNo tags found from file2\n", RESET;
    }
    else {
        print BOLD CYAN, "\n\nList of tags from file2: " . "@tagsArr2\n", RESET;
    }

    # simple intersection and union without List::Compare
    my %seen1 = map { $_ => 1 } @tagsArr;
    my %seen2 = map { $_ => 1 } @tagsArr2;

    my %union = (%seen1, %seen2);
    my @union = sort keys %union;

    my @intersection;
    for my $tag (@union) {
        push @intersection, $tag if $seen1{$tag} && $seen2{$tag};
    }

    if ( scalar @intersection == 0 ) {
        print BOLD CYAN, "\nNo intersection found between files\n", RESET;
    }
    else {
        print BOLD CYAN,
          "\nIntersection between file1 and file2: " . "@intersection\n", RESET;
    }

    if ( scalar @union != 0 ) {
        print BOLD CYAN, "\nUnion between file1 and file2: " . "@union\n",
          RESET;
    }

    close $fh;
    close $fh2;
}

1;
