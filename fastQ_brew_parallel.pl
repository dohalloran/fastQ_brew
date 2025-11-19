#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use File::Temp qw(tempdir);
use POSIX ":sys_wait_h";

my $SCRIPT = 'fastQ_brew.pl';

my %opt = (
    jobs            => 4,         # number of parallel processes
    reads_per_chunk => 100_000,   # reads per chunk (4 lines per read)
    o               => 'filtered.fq',
    lib             => 'sanger',
    qf              => 0,
    lf              => 0,
    truseq          => 0,
    trim_l          => 0,
    trim_r          => 0,
    fasta           => 0,
    dna_rna         => 0,
    rev_comp        => 0,
    plex            => 0,
);

GetOptions(
    'i=s'        => \$opt{i},       # input FASTQ (required)
    'o=s'        => \$opt{o},       # final output file
    'lib=s'      => \$opt{lib},
    'qf=f'       => \$opt{qf},
    'lf=i'       => \$opt{lf},
    'truseq!'    => \$opt{truseq},
    'trim_l=i'   => \$opt{trim_l},
    'trim_r=i'   => \$opt{trim_r},
    'fasta!'     => \$opt{fasta},
    'dna_rna!'   => \$opt{dna_rna},
    'rev_comp!'  => \$opt{rev_comp},
    'plex!'      => \$opt{plex},    # NOTE: ignored for parallel chunks
    'jobs=i'     => \$opt{jobs},
    'reads_per_chunk=i' => \$opt{reads_per_chunk},
) or die "Error in command line arguments\n";

die "Input FASTQ required (--i)\n" unless defined $opt{i};

my $input  = $opt{i};
my $output = $opt{o};

#----------------------#
#    Split into chunks #
#----------------------#

my $tmpdir = tempdir( "fastq_brew_parallelXXXXXX", CLEANUP => 1 );
my @chunks;

{
    open my $fh, '<', $input or die "Cannot open $input: $!";

    my $chunk_idx      = 0;
    my $reads_in_chunk = 0;
    my $chunk_fh;
    my $lines_per_read = 4;
    my $reads_per_chunk = $opt{reads_per_chunk};

    while (1) {
        my @lines;
        for (my $i = 0; $i < $lines_per_read; $i++) {
            my $line = <$fh>;
            last unless defined $line;
            push @lines, $line;
        }
        last unless @lines == $lines_per_read; # EOF

        if (!$chunk_fh || $reads_in_chunk >= $reads_per_chunk) {
            close $chunk_fh if $chunk_fh;
            $chunk_idx++;
            $reads_in_chunk = 0;
            my $chunk_file = "$tmpdir/chunk_$chunk_idx.fastq";
            open $chunk_fh, '>', $chunk_file or die "Cannot open $chunk_file: $!";
            push @chunks, $chunk_file;
        }

        print $chunk_fh @lines;
        $reads_in_chunk++;
    }

    close $chunk_fh if $chunk_fh;
    close $fh;
}

die "No reads found in input file $input\n" unless @chunks;

print "Created ", scalar(@chunks), " chunks in $tmpdir\n";

#----------------------#
#   Run in parallel    #
#----------------------#

my @pids;
my $max_jobs = $opt{jobs};

for my $chunk (@chunks) {
    while (@pids >= $max_jobs) {
        my $done = waitpid(-1, 0);
        @pids = grep { $_ != $done } @pids;
    }

    my $pid = fork();
    die "fork() failed: $!" unless defined $pid;

    if ($pid == 0) {
        # child
        my $chunk_out = "$chunk.out";

        my @cmd = (
            $^X, $SCRIPT,
            '--i', $chunk,
            '--o', $chunk_out,
            '--lib', $opt{lib},
            '--qf', $opt{qf},
            '--lf', $opt{lf},
            '--trim_l', $opt{trim_l},
            '--trim_r', $opt{trim_r},
        );

        push @cmd, '--truseq'   if $opt{truseq};
        push @cmd, '--fasta'    if $opt{fasta};
        push @cmd, '--dna_rna'  if $opt{dna_rna};
        push @cmd, '--rev_comp' if $opt{rev_comp};

        # NOTE: we do NOT forward --plex here (it compares two full files)

        print "Running: @cmd\n";
        exec @cmd or die "exec(@cmd) failed: $!";
    }
    else {
        push @pids, $pid;
    }
}

# wait for remaining children
while (@pids) {
    my $done = waitpid(-1, 0);
    @pids = grep { $_ != $done } @pids;
}

#----------------------#
#  Concatenate output  #
#----------------------#

open my $out_fh, '>', $output or die "Cannot open final output $output: $!";

for my $chunk (@chunks) {
    my $chunk_out = "$chunk.out";
    open my $cfh, '<', $chunk_out or die "Cannot open $chunk_out: $!";
    while (my $line = <$cfh>) {
        print $out_fh $line;
    }
    close $cfh;
}

close $out_fh;

print "Parallel fastQ_brew finished. Output written to $output\n";
