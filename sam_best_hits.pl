#' Author: Youtao Lu<luyoutao@sas.upenn.edu>
#' Copyright (c) 2020 Kim Laboratory, University of Pennsylvania
#' All Rights Reserved

use strict;
use warnings;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;

our $VERSION = 0.71;
our $LOGGER  = get_logger(__PACKAGE__);
my ( $inFiles, $labels ) = ( "", "" );
my @inFiles;
my @labels;
my ( $outFile, $outFh ) = ( "", undef );
my $endType     = "PE";
my $nameSorted  = 0;
my $scheme      = "AS";
my $as_per_mate = 0;
my $nm_per_mate = 0;
my $penalty     = 2;
my $ncores      = 1;
my $help        = 0;
my $version     = 0;
my $debug       = "info";

sub usage {
    print STDERR <<DOC;
Summary:
    Compare primary alignment scores from several SAM/BAM files. For each read (read pair if paired-end [PE]), output score from the respective input (blank if absent) as well as which one has the best score (multiple if ties). Three ways of calculating the score: 
    *) AS: original score reported in the tag 'AS';
    *) heuristic_raw: raw mapped length - penalty * #. of mismatches;
    *) heuristic_nonoverlap: nonoverlapping mapped length - penalty * #. of mismatches.

Usage:
    perl $0 --inFiles <input1.bam,input2.bam,...> --labels <species1,species2,...> [--outfile output.txt]  [--endType PE] [--nameSorted] [--ncores 1] [--scheme AS] [--AS-per-mate] [--NM-per-mate] [--penalty 2] [--debug info] [--version]

Options:
    --inFiles, -i   input files, comma separated, can be SAM or BAM format;
    --labels        labels of input files, comma separated;
    --outFile, -o   output file; if omitted, write to STDOUT; otherwise, if ending with '.gz', will be GZ compressed;
    --endType       choose from SE or PE (default: PE);
    --nameSorted    whether the input files are ALL name sorted. If not, they will be name sorted and saved to
                    temporary files with '.nameSorted.bam' suffix;
    --ncores        number of cores to use for name sorting (default: 1);
    --scheme        score scheme, choose from AS, heuristic_raw and heuristic_nonoverlap (default: AS);
    --AS-per-mate   boolean switch, if on, a read pair's AS equals to the sum of two mates' AS. This should be 
                    turned off if STAR, otherwise if Bowtie2 or BWA inputs (default: off); 
    --NM-per-mate   boolean switch, if on, instead of 'nM', a read pair's #. of mismatches equals to the sum of 
                    two mates' 'NM' (default: off)--this requires presence of tag 'NM'; note that STAR by default
                    does NOT output 'NM' (default: off); 
    --penalty       linear panelty of # of mismatches in heuristic scores (default: 2);
    --debug         level of logging, choose from fatal, error, warn, info, debug, trace (default: info);
    --version, -v   output version number.

Output:
    The output has n+3 columns where n is the number of input files/labels:
        *) read ID;
        *) score from 1st input;
        ...
        *) score from n-th input;
        *) best score;
        *) source label(s) of the best score(s);

Example:
    perl sam_best_hits.pl --inFiles human.sam,mouse.sam,bacteria.sam --labels human,mouse,bacteria 2>err
    ReadID  human   mouse   bacteria        BestScore       BestLabel
    NB501328:197:HMK3KBGX7:1:11101:1498:13562       49      41              49      human
    NB501328:197:HMK3KBGX7:1:11101:1744:6598        39      47              47      mouse
    NB501328:197:HMK3KBGX7:1:11101:2468:19747       105     44              105     human
    NB501328:197:HMK3KBGX7:1:11101:2835:13039       35      35      31      35      human,mouse
    NB501328:197:HMK3KBGX7:1:11101:2927:13134               78      44      78      mouse
    NB501328:197:HMK3KBGX7:1:11101:3669:18422               105     107     107     bacteria
DOC
}

GetOptions(
    "inFiles|i=s" => \$inFiles,
    "labels=s"    => \$labels,
    "outFile|o=s" => \$outFile,
    "endType=s"   => \$endType,
    "nameSorted"  => \$nameSorted,
    "ncores=i"    => \$ncores,
    "scheme=s"    => \$scheme,
    "AS-per-mate" => \$as_per_mate,
    "NM-per-mate" => \$nm_per_mate,
    "penalty=i"   => \$penalty,
    "help|h"      => \$help,
    "version|v"   => \$version,
    "debug=s"     => \$debug,
) or print STDERR "Wrong arguments!\n" && &usage() && exit(-1);

print "$0 v$VERSION\n" && exit(0) if $version;
&usage() && exit(0) if $help;

die("--inFiles empty!\n") if !defined($inFiles);
die("--labels empty!\n")  if !defined($labels);
@inFiles = split( ",", $inFiles );
@labels  = split( ",", $labels );
die("--inFiles do not match --labels!\n") if $#inFiles != $#labels;
for my $f ( grep { !-e $_ || $_ !~ /\.(sam|bam)$/i } (@inFiles) ) {
    die("$f does not exist or is not a SAM|BAM!\n");
}
my %seen;
for my $l (@labels) {
    next if !$seen{$l}++;
    die("$l is found duplicate in --labels!\n");
}
die("--endType can only be one of 'PE' and 'SE'!\n")
  if $endType ne "PE" && $endType ne "SE";
die(
"--scheme can only be one of 'AS', 'heuristic_raw' and 'heuristic_nonoverlap'!\n"
  )
  if $scheme ne "AS"
  && $scheme ne "heuristic_raw"
  && $scheme ne "heuristic_nonoverlap";
die("--penalty cannot be negative!\n") if $penalty < 0;

if ( $debug eq "fatal" ) {
    $LOGGER->level($FATAL);
}
elsif ( $debug eq "error" ) {
    $LOGGER->level($ERROR);
}
elsif ( $debug eq "warn" ) {
    $LOGGER->level($WARN);
}
elsif ( $debug eq "info" ) {
    $LOGGER->level($INFO);
}
elsif ( $debug eq "debug" ) {
    $LOGGER->level($DEBUG);
}
elsif ( $debug eq "trace" ) {
    $LOGGER->level($TRACE);
}
else {
    die("--debug can only be one of fatal, error, warn, info, debug, trace!\n");
}

my $appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::Screen");
my $layout   = Log::Log4perl::Layout::PatternLayout->new(
    "[%d{yyyy-MM-dd HH:mm:ss.SSS Z} %p] %m");
$appender->layout($layout);
$LOGGER->add_appender($appender);

$LOGGER->info(
"{ VERSION = $VERSION, inFiles = $inFiles, labels = $labels, outFile = $outFile, endType == $endType, nameSorted = $nameSorted, ncores = $ncores, scheme = $scheme, AS-per-mate = $as_per_mate, NM-per-mate = $nm_per_mate, penalty = $penalty, help = $help, debug = $debug, version = $version }\n"
);

{

    package Read;
    my $endType;

    sub import {
        my $class = shift;
        $endType = shift;
    }

    sub new {
        my $class = shift;
        my ( $readID, $flag, $chr, $pos, $cigar, @tags ) = @_;

        my $self = bless {
            readID       => $readID,
            flag         => $flag,
            chr          => $chr,
            pos          => $pos,
            whichInPair  => undef,
            isRev        => undef,
            cigar        => $cigar,
            tags         => \@tags,
            mappedLength => undef,
            NM           => undef,
            nM           => undef,
            AS           => undef,
        }, $class;
        return $self;
    }

    sub parse_flag {
        my $self = shift;
        $self->{whichInPair} =
          $endType eq "PE" ? ( $self->{flag} & 0x40 ? "R1" : "R2" ) : "R1";
        $self->{isRev} = $self->{flag} & 0x10 ? 1 : 0;
    }

    sub parse_cigar {
        my $self         = shift;
        my $mappedLength = 0;
        while ( $self->{cigar} =~ s/(\d+)[MX=]// ) {
            $mappedLength += $1;
        }
        $self->{mappedLength} = $mappedLength;
    }

    sub parse_NM {
        my $self = shift;
        my $NM;
        my @t = grep { /^NM/ } @{ $self->{tags} };
        $NM = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{NM} = $NM;
    }

    sub parse_nM {
        my $self = shift;
        my $nM;
        my @t = grep { /^nM/ } @{ $self->{tags} };
        $nM = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{nM} = $nM;
    }

    sub parse_AS {
        my $self = shift;
        my $AS;
        my @t = grep { /^AS/ } @{ $self->{tags} };
        $AS = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{AS} = $AS;
    }

}

{

    package ReadPair;
    my ( $endType, $scheme, $penalty, $as_per_mate, $nm_per_mate );

    sub import {
        my $class = shift;
        ( $endType, $scheme, $penalty, $as_per_mate, $nm_per_mate ) = @_;
    }

    sub new {
        my $class = shift;
        my $pair  = shift;
        my $self  = bless {
            R1               => $pair->{R1},
            R2               => $pair->{R2},
            readID           => undef,
            rawLength        => undef,
            nonoverlapLength => undef,
            mismatches       => undef,
            score            => undef,
        }, $class;
        $self->{readID} =
          defined( $self->{R1} )
          ? $self->{R1}->{readID}
          : $self->{R2}->{readID};
        return $self;
    }

    sub is_discordant {
        my $self = shift;
        $self->{R1}->{chr} ne $self->{R2}->{chr} ? 1 : 0;
    }

    sub parse_mappedLength {
        my $self = shift;
        if ( $endType eq "PE" ) {
            if ( defined( $self->{R1} ) && defined( $self->{R2} ) ) {
                $self->{R1}->parse_cigar();
                $self->{R2}->parse_cigar();
                $self->{rawLength} =
                  $self->{R1}->{mappedLength} + $self->{R2}->{mappedLength};
                if ( $self->is_discordant() ) {
                    $self->{nonoverlapLength} = $self->{rawLength};
                    return;
                }
                my $overlap =
                    $self->{R1}->{isRev}
                  ? $self->{R2}->{pos} + $self->{R2}->{mappedLength} -
                  $self->{R1}->{pos}
                  : $self->{R1}->{pos} +
                  $self->{R1}->{mappedLength} -
                  $self->{R2}->{pos};
                $overlap = $overlap > 0 ? $overlap : 0;
                $self->{nonoverlapLength} = $self->{rawLength} - $overlap;
            }
            elsif ( defined( $self->{R1} ) ) {
                $self->{R1}->parse_cigar();
                $self->{rawLength} = $self->{nonoverlapLength} =
                  $self->{R1}->{mappedLength};
            }
            else {
                $self->{R2}->parse_cigar();
                $self->{rawLength} = $self->{nonoverlapLength} =
                  $self->{R2}->{mappedLength};
            }
        }
        else {
            $self->{R1}->parse_cigar();
            $self->{rawLength} = $self->{nonoverlapLength} =
              $self->{R1}->{mappedLength};
        }
    }

    sub parse_mismatches {
        my $self = shift;
        if ( $endType eq "PE" ) {
            if ($nm_per_mate) {
                my $s1 = defined( $self->{R1} ) ? $self->{R1}->parse_NM() : 0;
                my $s2 = defined( $self->{R2} ) ? $self->{R2}->parse_NM() : 0;
                $self->{mismatches} = $s1 + $s2;
            }
            else {
                $self->{mismatches} =
                  defined( $self->{R1} )
                  ? $self->{R1}->parse_nM()
                  : $self->{R2}->parse_nM();
            }
        }
        else {
            if ($nm_per_mate) {
                $self->{mismatches} = $self->{R1}->parse_NM();
            }
            else {
                $self->{mismatches} = $self->{R1}->parse_nM();
            }
        }
    }

    sub parse_AS {
        my $self = shift;
        my @t;
        if ( $endType eq "PE" ) {
            if ($as_per_mate) {
                my $s1 = defined( $self->{R1} ) ? $self->{R1}->parse_AS() : 0;
                my $s2 = defined( $self->{R2} ) ? $self->{R2}->parse_AS() : 0;
                $self->{score} = $s1 + $s2;
            }
            else {
                $self->{score} =
                  defined( $self->{R1} )
                  ? $self->{R1}->parse_AS()
                  : $self->{R2}->parse_AS();
            }
        }
        else {
            $self->{score} = $self->{R1}->parse_AS();
        }
    }

    sub cal_score {
        my $self = shift;
        if ( $scheme eq "AS" ) {
            $self->parse_AS();
        }
        else {
            $self->parse_mappedLength();
            $self->parse_mismatches();
            if ( $scheme eq "heuristic_raw" ) {
                $self->{score} =
                  $self->{rawLength} - $penalty * $self->{mismatches};
            }
            else {
                $self->{score} =
                  $self->{nonoverlapLength} - $penalty * $self->{mismatches};
            }
        }
    }
}

{

    package SamReader;

    sub new {
        my $class = shift;
        my ( $inFile, $label, $endType, $nameSorted, $ncores ) = @_;
        my $self = bless {
            inFile     => $inFile,
            inFh       => undef,
            label      => $label,
            endType    => $endType,
            nameSorted => $nameSorted,
            ncores     => $ncores,
            read       => undef,
            readPair   => undef,
        }, $class;
        return $self;
    }

    sub init_fh {
        my $self = shift;

        if ( !$self->{nameSorted} ) {
            my $tmpFile = $self->{inFile};
            $tmpFile =~ s/\.(bam|sam)$//i;
            $tmpFile .= ".nameSorted" . ".bam";
            $LOGGER->warn(
"$self->{inFile} is not name sorted. Name sorting it and saving to $tmpFile...\n"
            );
            $LOGGER->warn("$tmpFile exists already! It will be overwritten.\n")
              if -e $tmpFile;
            my $cmd =
"samtools sort -n -\@ $self->{ncores} -o $tmpFile $self->{inFile}";
            $LOGGER->debug("\$cmd = $cmd\n");
            my $exit = system($cmd);
            $LOGGER->fatal("Error happended when sorting $self->{inFile}!\n")
              && die()
              if $exit;
            $self->{inFile} = $tmpFile;
        }
        $LOGGER->info("Opening file handle for $self->{inFile}...\n");
        open( $self->{inFh}, "samtools view -F 0x100 $self->{inFile} |" )
          or $LOGGER->fatal("Cannot open $self->{inFile} for read!\n") && die();
    }

    sub fin_fh {
        my $self = shift;
        $LOGGER->info("Closing file handle for $self->{inFile}...\n");
        close $self->{inFh} or $LOGGER->warn("Cannot close $self->{inFile}!\n");
    }

    sub get_read {
        my $self = shift;
        my $inFh = $self->{inFh};
        my $line = <$inFh>;
        return undef if !defined($line);
        chomp($line);
        my ( $readID, $flag, $chr, $pos, $cigar, @tags, @F, $read );
        @F      = split( "\t", $line );
        $readID = $F[0];
        $flag   = $F[1];
        $chr    = $F[2];
        $pos    = $F[3];
        $cigar  = $F[5];
        @tags   = @F[ 11 .. $#F ];
        $read   = Read->new( $readID, $flag, $chr, $pos, $cigar, @tags );
        $read->parse_flag();
        return $read;
    }

    # @Get a read (pair), compute its score
    sub next {
        my $self = shift;
        my ( $read, $readPair );
        $read = $self->get_read();
        if ( $self->{endType} eq "PE" ) {
            if ( !defined($read) ) {    # has reached end of file
                if ( defined( $self->{read} ) ) {  # we have a residual in cache
                    $readPair = ReadPair->new(
                        { $self->{read}->{whichInPair} => $self->{read} } );
                    $self->{read} = undef;
                    $readPair->cal_score();
                    $self->{readPair} = $readPair;
                    $LOGGER->trace(
"A \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                    );
                }
                else {                             # empty cache
                    $self->{readPair} = undef;
                }
            }
            else {
                if ( defined( $self->{read} ) ) {    # the cache is not empty
                    if ( $self->{read}->{readID} eq $read->{readID} )
                    {    # current read pairs with last one
                        $readPair = ReadPair->new(
                            {
                                $self->{read}->{whichInPair} => $self->{read},
                                $read->{whichInPair}         => $read,
                            }
                        );
                        $readPair->cal_score();
                        $self->{readPair} = $readPair;
                        $LOGGER->trace(
"C \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                        );
                        $self->{read} = undef;
                    }
                    else
                    { # current read is different from last one, process last one, cache the current one
                        $readPair = ReadPair->new(
                            { $self->{read}->{whichInPair} => $self->{read} } );
                        $readPair->cal_score();
                        $self->{readPair} = $readPair;
                        $LOGGER->trace(
"D \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                        );
                        $self->{read} = $read;
                    }
                }
                else
                { # empty cache: last was a read pair and has been processed, cache the current read, then get a new read
                    $self->{read} = $read;
                    $read = $self->get_read();
                    if ( !defined($read) ) {    # has reached end of file
                        $readPair = ReadPair->new(
                            { $self->{read}->{whichInPair} => $self->{read} } );
                        $readPair->cal_score();
                        $self->{readPair} = $readPair;
                        $LOGGER->trace(
"E \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                        );
                        $self->{read} = undef;
                    }
                    else {                      # new read exists
                        if ( $self->{read}->{readID} eq $read->{readID} )
                        {    # current read pairs with last one
                            $readPair = ReadPair->new(
                                {
                                    $self->{read}->{whichInPair} =>
                                      $self->{read},
                                    $read->{whichInPair} => $read,
                                }
                            );
                            $readPair->cal_score();
                            $self->{readPair} = $readPair;
                            $LOGGER->trace(
"F \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                            );
                            $self->{read} = undef;
                        }
                        else
                        { # current read is different from last one, process the cached read as a singleton, then cache the current read
                            $readPair = ReadPair->new(
                                {
                                    $self->{read}->{whichInPair} =>
                                      $self->{read}
                                }
                            );
                            $readPair->cal_score();
                            $self->{readPair} = $readPair;
                            $LOGGER->trace(
"G \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                            );
                            $self->{read} = $read;
                        }
                    }
                }
            }
        }
        else {    # SE
            if ( !defined($read) ) {    # has reached end of file
                $self->{readPair} = undef;
            }
            else {                      # not end of file yet
                $readPair = ReadPair->new( { R1 => $read } );
                $readPair->cal_score();
                $self->{readPair} = $readPair;
                $LOGGER->trace(
"H \$self->{readPair}->{readID} = $self->{readPair}->{readID}\n"
                );
            }
        }
    }
}

{

    package SamReaders;
    use IO::Zlib;
    use List::Util qw( reduce max );
    use Inline 'C' => <<CODE;
/* from samtools/1.9/bam_sort.c */
int strnum_cmp(const char *_a, const char *_b)
{
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}
CODE

    sub new {
        my $class = shift;
        my (
            $inFiles,     $labels, $outFile, $endType,
            $nameSorted,  $ncores, $scheme,  $penalty,
            $as_per_mate, $nm_per_mate
        ) = @_;
        Read->import($endType);    # Initialize class static
        ReadPair->import( $endType, $scheme, $penalty, $as_per_mate,
            $nm_per_mate );        # Initialize class static
        my $self = bless {
            inFiles    => $inFiles,
            labels     => $labels,
            samReaders => [],
            outFile    => $outFile,
            outFh      => undef,
            endType    => $endType,
            nameSorted => $nameSorted,
            ncores     => $ncores,
            penalty    => $penalty,
        }, $class;
        return $self;
    }

    sub init {
        my $self = shift;
        my $n    = scalar @{ $self->{inFiles} };
        $self->{samReaders} = [
            map {
                SamReader->new(
                    ${ $self->{inFiles} }[$_],
                    ${ $self->{labels} }[$_],
                    $self->{endType}, $self->{nameSorted}, $self->{ncores}
                )
            } ( 0 .. ( $n - 1 ) )
        ];
        if ( $self->{outFile} eq '' ) {
            $LOGGER->info("Opening file handle for STDOUT...\n");
            $self->{outFh} = *STDOUT;
        }
        elsif ( $self->{outFile} =~ /\.(gz|gzip)$/i ) {
            $LOGGER->info("Opening file handle for $self->{outFile}...\n");
            ( $self->{outFh} = IO::Zlib->new( $self->{outFile}, 'w' ) )
              or $LOGGER->fatal("Cannot open $self->{outFile} for write!\n")
              && die();
        }
        else {
            $LOGGER->info("Opening file handle for $self->{outFile}...\n");
            open( $self->{outFh}, ">", $self->{outFile} )
              or $LOGGER->fatal("Cannot open $self->{outFile} for write!\n")
              && die();
        }

        $LOGGER->debug( "({SamReaders->init}) \@{\$self->{samReaders}} = ",
            scalar @{ $self->{samReaders} }, "\n" );

        for my $samReader ( @{ $self->{samReaders} } ) {
            $samReader->init_fh();
            $samReader->next();    # get the first readPair
            $LOGGER->trace(
"\$samReader->{readPair}->{readID} = " . (defined($samReader->{readPair}) ? $samReader->{readPair}->{readID} : "") . "\n"
            );
        }
        my $outFh  = $self->{outFh};
        my @labels = @{ $self->{labels} };
        print $outFh
          join( "\t", ( "ReadID", @labels, "BestScore", "BestLabel" ) ), "\n";
    }

    sub iter {
        my $self   = shift;
        my $n      = scalar @{ $self->{inFiles} };
        my $labels = $self->{labels};
        my $outFh  = $self->{outFh};
        while (1) {
            my $samReaders = $self->{samReaders};
            my @idcs_hasRead =
              grep { defined( ${$samReaders}[$_]->{readPair} ) }
              ( 0 .. ( $n - 1 ) );
            $LOGGER->debug("\@idcs_hasRead = @idcs_hasRead\n");
            last if !scalar @idcs_hasRead;
            my @readIDs =
              map { ${$samReaders}[$_]->{readPair}->{readID} } @idcs_hasRead;
            $LOGGER->debug("\@readIDs = @readIDs, \$#readIDs = $#readIDs\n");
            $LOGGER->debug("\@idcs_hasRead = @idcs_hasRead\n");
            my $minReadID =
              @readIDs > 1
              ? reduce { &strnum_cmp( $a, $b ) < 0 ? $a : $b } @readIDs
              : shift @readIDs;
            $LOGGER->debug("\$minReadID = $minReadID\n");
            my @idcs_minReadID =
              grep { ${$samReaders}[$_]->{readPair}->{readID} eq $minReadID }
              @idcs_hasRead;
            $LOGGER->debug("\@idcs_minReadID = @idcs_minReadID\n");
            my @labels_minReadID = @{$labels}[@idcs_minReadID];
            $LOGGER->debug("\@labels_minReadID = @labels_minReadID\n");
            my @scores =
              map { ${$samReaders}[$_]->{readPair}->{score} } (@idcs_minReadID);
            $LOGGER->debug("\@scores = @scores\n");
            my $maxScore = max(@scores);
            $LOGGER->debug("\$maxScore = $maxScore\n");
            my @idcs_maxScore =
              grep { ${$samReaders}[$_]->{readPair}->{score} == $maxScore }
              (@idcs_minReadID);
            $LOGGER->debug("\@idcs_maxScore = @idcs_maxScore\n");
            my @labels_maxScore = @{$labels}[@idcs_maxScore];
            $LOGGER->debug("\@labels_maxScore = @labels_maxScore\n");
            my %lookup = map { ( $_, 1 ) } @idcs_minReadID;
            print $outFh ($minReadID);
            map {
                print $outFh (
                    "\t",
                    defined( $lookup{$_} )
                    ? ${$samReaders}[$_]->{readPair}->{score}
                    : ""
                );
            } ( 0 .. ( $n - 1 ) );
            print $outFh "\t", $maxScore;
            print $outFh "\t", join( ",", @labels_maxScore );
            print $outFh "\n";

            for my $i (@idcs_minReadID) {
                $LOGGER->debug("Updating readPair for $labels[$i]\n");
                ${$samReaders}[$i]->next();
            }
        }
    }

    sub fin {
        my $self = shift;
        $LOGGER->info("Closing output file handle...\n");
        close $self->{outFh}
          or
          $LOGGER->warn("Failed to close file handle for $self->{outFile}!\n");

        for my $samReader ( @{ $self->{samReaders} } ) {
            $samReader->fin_fh();
        }
    }
}

my $samReaders = SamReaders->new(
    \@inFiles, \@labels, $outFile, $endType,     $nameSorted,
    $ncores,   $scheme,  $penalty, $as_per_mate, $nm_per_mate
);
$samReaders->init();
$samReaders->iter();
$samReaders->fin();
$LOGGER->info("All done!\n");
