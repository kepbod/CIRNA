#!/usr/bin/perl
use strict;
use warnings;
use Carp;
use File::Basename;
use List::Util qw(sum);
use lib '/picb/rnomics1/xiaoou/program/usefullib/perl/';
use Overlap qw(OverlapMax OverlapMap OverlapMerge);
use SAMParse qw(ReadSplit);
use RefFileParse qw(ExtractInfo);

#
# Program Name: ILRFind.pl
# Function: Find Intron-derived LncRNA seed regions.
#

our $AUTHOR = "Xiao-Ou Zhang";
our $VERSION = "1.5.0";

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

# set essential messages
sub HELP_MESSAGE {
    print <<HELP;
Usage: ILRFind.pl -r <reference file> -b <BAM file> -s <reads number> [-i] [-q]
                 [-w <window cutoff>] [-e <extend cutoff>] [-a <added ratio cutoff>]
                 [-l <read length>] [-n <chromosome size>] [-c <chromosome>]
                 [-d <output directory>] [-o <output file>]
Notice: 1 The BAM file must contain '.bam' and the letters before them will
          be used as the name of output, log and error file.
        2 The default values of each optional parameter:
          window cutoff(5), extend cutoff(3), added ratio cutoff(0.3),
          read lenth(75), chromosome size(24)
        3 The criterion used in this program:
          (1) length of each region is longer than 200bp
          (2) length of starting window is 100bp
          (3) number of bins in each starting window reaching the window cutoff
              (default if 5) is more than 4 (there are 5 bins in each starting window)
          (4) when extending window, the BPKM of added bin is larger than 30% of
              average BPMK of the window, and larger than extend cutoff (default is 3)
        4 If the 'i' is set, the maximum intron setting will be used.
        5 If the 'q' is set, the statistics mode will be set.
        6 The program will creat a log file and a error file. The log file will
          record the running info, and the error file will record error info. If
          there are no errors when running this program, the error file will be
          deleted automatically.
HELP
}
sub VERSION_MESSAGE {
    print <<VERSION;
ILRFind - Find Intron-derived LncRNA seed regions.
Version: $VERSION, Maintainer: $AUTHOR
VERSION
}
sub help {
    print <<SIMPLE_HELP;
Usage: ILRFind.pl -r <reference file> -b <BAM file> -s <reads number> [-i] [-q]
                 [-w <window cutoff>] [-e <extend cutoff>] [-a <added ratio cutoff>]
                 [-l <read length>] [-n <chromosome size>] [-c <chromosome>]
                 [-d <output directory>] [-o <output file>]
See 'ILRFind.pl --help' for more details.
SIMPLE_HELP
    exit;
}

# get options and check
my %opts = ( w=>'5', e=>'3', a=>'0.3', l=>'75', n=>'24' );
getopts('hiqr:b:s:w:e:a:l:n:c:d:o:', \%opts) or help();
help() if $opts{h} or not exist_essential_parameter();

# check and set output directory
if(exists $opts{d}) {
    $opts{d} =~ s/([^\/]$)/$1\//;
    -d $opts{d} or croak "The diresctor you input is wrong!";
}
else {
    $opts{d} = '';
}

# check and set output filename
my $bamname = fileparse($opts{b});
$bamname =~ /(.*)\.bam/;
$opts{o} = $1 if not exists $opts{o};

# set output, log and error file
open my $result,'>',"$opts{d}$opts{o}.id-LncRNA" or croak "Can't open file $opts{d}$opts{o}!";
open STDOUT,'>',"$opts{d}$opts{o}.log";
open STDERR,'>',"$opts{d}$opts{o}.error";
my $st;
if ($opts{q}) {
    open $st,'>',"$opts{d}$opts{o}.statistics" or croak "Can't open file $opts{d}$opts{o}!"
}

# set cutoffs
my $window_cutoff = $opts{w} * $opts{l} * $opts{s} * 20 / 10**9;
my $extend_cutoff = $opts{e} * $opts{l} * $opts{s} * 20 / 10**9;

# set chromosome
my $escape = '';
if ($opts{n} =~ /([XY])$/) {
    $escape = $1;
    $opts{n} =~ s/[XY]$//;
    ++$opts{n};
}
my @chrom;
if (not exists $opts{c}) {
    for (1..$opts{n}-2) {
        push @chrom,'chr' . $_;
    }
    if ($escape eq 'X') {
        push @chrom,'chrY'
    }
    elsif ($escape eq 'Y') {
        push @chrom,'chrX'
    }
    else {
        push @chrom,'chrX';
        push @chrom,'chrY';
    }
}
else {
    @chrom = ($opts{c});
}

# record command and starting time
my $command;
while (my ($key, $value) = each %opts) {
    $command .= " -$key $value";
}
print "Command: $0$command\n";
print 'Program starts at ', scalar(localtime(time)), "\n";

# write table header
print $result "Position\tLength\tBPKM\tGene\n";

for my $chr (@chrom) { # for each chromosome
    # record running time
    print "Start $chr at time ", scalar(localtime(time)), "\n";

    # calculate intron regions
    my $intron_index_ref;
    {
        # read reference file
        open my $ref_file,'<',$opts{r} or croak "Can't open $opts{r} file!";
        my $refname = fileparse($opts{r}, qr/\Q.txt\E/);
        my ($gene_sta_end_ref, $exon_sta_end_ref, $intron_sta_end_ref)
            = (ExtractInfo($ref_file, $refname, $chr, '-', 1))[0,2,3];
        close $ref_file;

        if ($opts{i}) {
            # calculate intron regions
            $intron_index_ref = OverlapMax($intron_sta_end_ref, '-', 1);
        }
        else {
            # calculate gene regions
            my $max_gene_ref = OverlapMax($gene_sta_end_ref, '-', 1);

            # calculate exon regions
            my $max_exon_ref = OverlapMax($exon_sta_end_ref, '-');

            # calculate extra_exon regions (contain intron and non-gene regions)
            my @extra_exon;
            my $sta = (split /-/,$max_exon_ref->[0])[1]; # the higher boundary
            my $end; # the lower boundary
            for (1..$#{$max_exon_ref}) {
                $end = (split /-/,$max_exon_ref->[$_])[0];
                push @extra_exon,$sta . '-' . $end;
                $sta = (split /-/,$max_exon_ref->[$_])[1];
            }

            # merge gene regions and extra_exon regions to get intron regions
            $intron_index_ref = OverlapMerge($max_gene_ref,\@extra_exon,'-', 10);
        }
    }

    LABEL:for my $intron (@$intron_index_ref) { # for each intron

        # length is not enough
        next if (split /-/,$intron)[1] - (split /-/,$intron)[0] < 200;

        # extract gene info
        my %gene;
        map {$gene{$_} = 1} split /-/,(split /-/,$intron,3)[-1];
        my $name = join " ",map { (my $x = $_) =~ s/#/-/g; $x } keys %gene;
        $intron = join "-",(split /-/,$intron)[0,1];

        # split intron into bins (20bp each)
        my $bins_ref = intron_split($intron);

        # calculate bases in each bin and flag to indicate whether reaching cutoff
        my ($base_in_bins_ref, $flag_ref);
        {
            # read reads in this intron region
            open my $read_iterator,'-|',"samtools view $opts{b} $chr:$intron";

            # split reads
            my @read_segments;
            while (<$read_iterator>) { # for each read in this intron
                my ($pos, $cigar) = (split)[3, 5]; # get position and cigar
                my $segments_ref = ReadSplit($pos, $cigar, '-'); # split read
                push @read_segments,@{$segments_ref};
            }
            close $read_iterator;

            next LABEL if not exists $read_segments[0]; # no reads

            # reads in each bin
            my $reads_map_to_bin_ref
                = OverlapMap($bins_ref, \@read_segments, '-', 10);

            # bases in each bin
            ($base_in_bins_ref, $flag_ref)
                = bin_to_base($reads_map_to_bin_ref, $#{$bins_ref}+1);
        }

        next if sum(@$flag_ref) < 4; # not enough bins reaching the cutoff

        my @slidewindows;
        my $window = 0;
        $window += $_ for @{$base_in_bins_ref}[0..4]; # first window
        if (sum(@{$flag_ref}[0..4]) >= 4) {
            # meet the criterion of starting window
            push @slidewindows,$window . ':' . 0;
        }
        for (5..$#{$base_in_bins_ref}) { # move slide window
            $window -= $base_in_bins_ref->[$_ - 5];
            $window += $base_in_bins_ref->[$_];
            if (sum(@{$flag_ref}[$_-4..$_]) >= 4) {
                # meet the criterion of starting window
                push @slidewindows,$window . ':' . ($_ - 4);
            }
        }

        next if not exists $slidewindows[0]; # no slide window

        # sort slide windows according to the base in each bin (descending order)
        my @candidate_windows = map { $_->[0] }
                                sort { $b->[1] <=> $a->[1] }
                                map { [ $_, (split /:/)[0] ] }
                                @slidewindows;

        my @candidate_LncRNA_windows;
        # initiate excluded bin set (start with all 0s)
        my $excluded_bins = '0' x ($#{$base_in_bins_ref} + 1);
        for (@candidate_windows) { # for each candidate window

            # make sure no windows belonging to the excluded bin set
            my $flag_of_excluded = substr $excluded_bins,(split /:/)[1],5;
            next if eval($flag_of_excluded);

            # extend window
            my ($extended_window, $bpkm, $flag, $statistics)
                = window_extend($_, $base_in_bins_ref, $excluded_bins, $opts{q});

            if ($opts{q}) {
                print $st $statistics;
            }

            if ($flag) { # if length of extended window is longer than 200bp
                push @candidate_LncRNA_windows,$extended_window . '-' . $bpkm;
                # refine excluded bin set (if excluded, set to 1)
                my ($left, $right) = split /-/,$extended_window;
                substr($excluded_bins, $left, ($right - $left + 1))
                    = '1' x ($right - $left + 1);
                }
            }

        # write result
        for (@candidate_LncRNA_windows) {
            my ($first, $last, $bpkm) = split /-/;
            my ($sta ,$end)
                = ((split /-/,$$bins_ref[$first])[0], (split /-/,$$bins_ref[$last])[1]);
            print $result $chr . ':' . $sta . '-' . $end . "\t"
                . ($end - $sta) . "\t" . $bpkm . "\t" . $name . "\n";
        }
    }
}

# close file handle
close $result;
if ($opts{q}) { close $st; }

# delete error file if no errors
unlink "$opts{d}$opts{o}.error" if -z "$opts{d}$opts{o}.error";

# record ending time
print 'Program ends at ', scalar(localtime(time)), "\n";

##########Subroutine##########

#
# Function: Check essential parameters.
#
sub exist_essential_parameter {
    my @essential_parameter = qw( r b s );
    exists $opts{$_} || return 0 for @essential_parameter;
    return 1;
}

#
# Function: Split intron into bins
#
sub intron_split {
    my ($pos1, $end) = split /-/,$_[0];
    my @bins;
    my $bin = 20;
    my $pos2 = $pos1 + $bin;
    my $i = 0;
    while ($pos2 <= $end) {
        push @bins,$pos1 . '-' . $pos2 . '-' . $i;
        $pos1 = $pos2;
        $pos2 += $bin;
        ++$i;
    }
    return \@bins;
}

#
# Function: Calculate base-pairs in each bin
#
sub bin_to_base {
    my $segments_ref = shift;
    my $bin_num = shift;
    my @base_in_bins = (0) x $bin_num;
    for (@{$segments_ref}) {
        my ($sta, $end, $index) = split /-/;
        $base_in_bins[$index] += $end - $sta;
    }
    my @flag = map {$base_in_bins[$_] >= $window_cutoff ? 1 : 0} 0..$bin_num-1;
    return (\@base_in_bins, \@flag);
}

#
# Function: extend slide window as much as possible
#
sub window_extend{
    my ($base, $n) = split /:/,$_[0];
    my $base_ref = $_[1];
    my $excluded_bins = $_[2];
    my $st_flag = $_[3];
    my $base_iterator = $base;
    my $bpkm = ($base * 10**9) / ($opts{l} * $opts{s} * 100);
    my $statistics = '';
    my ($i, $j) = (0, 0);
    while (1) {
        last if $n - $i <= 0;
        ++$i;
        if (substr($excluded_bins,$n - $i,1) == 1) {
            --$i;
            last;
        }
        if ($st_flag) {
            my $bin_bpkm = ($$base_ref[$n - $i] * 10**9) / ($opts{l} * $opts{s} * 20);
            $statistics .= "$bpkm\t$bin_bpkm\t$i\n";
        }
        $base_iterator += $$base_ref[$n - $i];
        my $new_bpkm
            = ($base_iterator * 10**9) / ($opts{l} * $opts{s} * (100 + 20 * $i));
        if ((not $st_flag and $new_bpkm < $bpkm * (5 + $i - 1 + $opts{a}) / (5 + $i))
                or $$base_ref[$n - $i] < $extend_cutoff) {
            $base_iterator -= $$base_ref[$n - $i];
            --$i;
            last;
        }
        $bpkm = $new_bpkm;
    }
    while (1) {
        last if $n + 4 + $j >= $#{$base_ref};
        ++$j;
        if (substr($excluded_bins,$n + 4 + $j,1) == 1) {
            --$j;
            last;
        }
        if ($st_flag) {
            my $bin_bpkm = ($$base_ref[$n + 4 + $j] * 10**9) / ($opts{l} * $opts{s} * 20);
            $statistics .= "$bpkm\t$bin_bpkm\t$j\n";
        }
        $base_iterator += $$base_ref[$n + 4 + $j];
        my $new_bpkm
            = ($base_iterator * 10**9) / ($opts{l} * $opts{s} * (100 + 20 * ($i + $j)));
        if ((not $st_flag and $new_bpkm < $bpkm * (5 + $i + $j - 1 + $opts{a}) / (5 + $i +$j))
                or $$base_ref[$n + 4 + $j] < $extend_cutoff) {
            $base_iterator -= $$base_ref[$n + 4 +$j];
            --$j;
            last;
        }
        $bpkm = $new_bpkm;
    }
    my $window = ($n - $i) . '-' . ($n + 4 + $j);
    my $flag = (5 + $i + $j) >= 10 ? 1 : 0;
    return ($window, $bpkm, $flag, $statistics);
}
