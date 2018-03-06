#
# File conversion module for fastQ_brew
#
# Please direct questions and support issues to <https://github.com/dohalloran/fastQ_brew/issues>
#
# Author: Damien O'Halloran, The George Washington University, 2017
#
# GNU GENERAL PUBLIC LICENSE
#
# POD documentation before the code

=head1 NAME

fastQ_brew_Conversions - file reformatting and demultiplex check for fastQ_brew

=head2 SYNOPSIS

  use Moose;
  use Modern::Perl;
  use base 'Exporter';

=head2 DESCRIPTION

This package provides subroutines for file conversion

=head2 Support

All contributions are welcome

=head2 Reporting Bugs

Report bugs to the fastQ_brew bug tracking system to help keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:
  https://github.com/dohalloran/fastQ_brew/issues

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package fastQ_brew_Conversions;

use Moose;
use namespace::autoclean;
use Modern::Perl;
use List::Compare;
use base 'Exporter';
use Term::ANSIColor qw(:constants);
use autodie;

our @EXPORT = qw/ _convert_fasta _rev_comp _dna_rna _de_plex /;

##################################

=head2 _convert_fasta()

 Title   : _convert_fasta()
 Usage   : _convert_fasta(%arg);
 Function: option to convert fastQ file to fastA
 Returns : fastA file
 Args    : %arg

=cut

##################################

sub _convert_fasta {
    my $in_file = shift;
    print BOLD CYAN, "\n\nConverting fastQ file to fastA...", RESET;
    my $temp;
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

        # Print to fasta file.
        print $fh_out ">$temp[0]\n";
        print $fh_out "$temp[1]\n";
    }
    close $fh;
    close $fh_out;
}
##################################

=head2 _reverse_comp()

 Title   : _reverse_comp()
 Usage   : $self->_reverse_comp(%arg)
 Function: option to rev comp fastQ reads
 Returns : reverse complemented fastQ file
 Args    : %arg

=cut

##################################

sub _rev_comp {
    my $in_file = shift;
    print BOLD CYAN, "\n\nReverse complementing fastQ reads...", RESET;
    my $temp;
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

        # rev comp the array element
        $temp[1] =~ tr/ATGCatgc/TACGtacg/;
        $revComp = reverse( $temp[1] );

        # Print to revcomp file.
        print $fh_out "$temp[0]\n";
        print $fh_out "$revComp\n";
        print $fh_out "$temp[2]\n";
        print $fh_out "$temp[3]\n";
    }
    close $fh;
    close $fh_out;
}
##################################

=head2 _dna_rna()

 Title   : _dna_rna()
 Usage   : $self->_dna_rna(%arg)
 Function: option to convert dna to rna for fastQ reads
 Returns : RNA fastQ file
 Args    : %arg

=cut

##################################

sub _dna_rna {
    my $in_file = shift;
    print BOLD CYAN, "\n\nTranscribing fastQ reads...", RESET;
    my $temp;
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

        # transcribe the array element
        $temp[1] =~ tr/Tt/Uu/;

        # Print to RNA file.
        print $fh_out "$temp[0]\n";
        print $fh_out "$temp[1]\n";
        print $fh_out "$temp[2]\n";
        print $fh_out "$temp[3]\n";
    }
    close $fh;
    close $fh_out;
}
####################################

=head2 _de_plex()

 Title   : _de_plex()
 Usage   : $self->_de_plex(%arg)
 Function: option to check that two FASTQ files were de-multiplexed correctly
 Returns : list of union and intersection tags from each file
 Args    : %arg

=cut

##################################

sub _de_plex {
    my (@fileArray) = @_;
    my $in_file  = $fileArray[0];
    my $in_file2 = $fileArray[1];
    print BOLD CYAN, "\n\nChecking fastQ files for proper demultiplexing...",
      RESET;
    my $temp;
    my @temp;

    open my $fh, '<', $in_file
      or die "Cannot open $in_file: $!";

    open my $fh2, '<', $in_file2
      or die "Cannot open $in_file2: $!";

    my @tagsArr;
    my @tagsArr2;
    while ( my $line = <$fh> ) {
        if ( $line =~ m/:([AGCT][AGCT][AGTC][AGCT][AGCT][AGTC])$/ ) {
            push( @tagsArr, $1 ) unless grep { $_ eq $1 } @tagsArr;
        }
    }
    while ( my $line2 = <$fh2> ) {
        if ( $line2 =~ m/:([AGCT][AGCT][AGTC][AGCT][AGCT][AGTC])$/ ) {
            push( @tagsArr2, $1 ) unless grep { $_ eq $1 } @tagsArr2;
        }
    }
    print BOLD CYAN, "\n\nList of tags from file1: ". "@tagsArr\n",  RESET;
    print BOLD CYAN, "\n\nList of tags from file2: ". "@tagsArr2\n", RESET;

    my $lc = List::Compare->new( \@tagsArr, \@tagsArr2 );

    my @intersection = $lc->get_intersection;
    my @union        = $lc->get_union;

    print BOLD CYAN, "\nIntersection between file1 and file2: " . "@intersection\n", RESET;
    print BOLD CYAN, "\nUnion between file1 and file2: " . "@union\n", RESET;

    close $fh;
    close $fh2;
}
####################################
####################################

1;
