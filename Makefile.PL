use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'fastQ_brew',
    AUTHOR           => q{Damien OHalloran <damienoh@gwu.edu>},
    VERSION_FROM     => 'lib/fastQ_brew.pm',
    ABSTRACT_FROM    => 'lib/fastQ_brew.pm',
    LICENSE          => 'gpl_2',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
        'Sys::Hostname'           => 0,
        'File::Copy'              => 0,
        'Term::ANSIColor'         => 0,
        'Moose'                   => 2.0,
        'Modern::Perl'            => 1.0,
        'Text::LevenshteinXS'     => 0,
        'List::MoreUtils'         => 0,
        'Win32::Console::ANSI'    => 0,
        'Statistics::Basic'       => 0,
        'List::Compare'           => 0
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'fastQ_brew-*' },
);
