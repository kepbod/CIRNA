#!/usr/bin/perl
use strict;
use warnings;
use Carp;
use File::Basename;
use List::Util qw(sum);
use lib '/picb/rnomics1/xiaoou/program/usefullib/perl/';
use Overlap qw(OverlapMap);
use SAMParse qw(ReadSplit);
use RefFileParse qw(ExtractInfo);

#
# Program Name: NEBPKM.pl
# Function: Calculate BPKM of id-LncRNA seed regions' neighbor exons.
#

our $AUTHOR = "Xiao-Ou Zhang";
our $VERSION = "1.2.0";

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

# set essential messages
sub HELP_MESSAGE {
    print <<HELP;
Usage: NEBPKM.pl -r <reference file> -b <BAM file> -f <id-LncRNA file> -s <read number>
                [-e] [-c] [-i] [-a additional BAM file] [-t additional BAM file read number]
                [-l <read length>] [-d <output directory>] [-o <output file>]
Notice: 1 If parameter 'e' was set, sole exon gene will be excluded.
        2 If parameter 'c' was set, compare mode will be used.
        3 If parameter 'a' was set, additional BPKM will be cauculated.
        4 If parameter 'i' was set, interesting flag will be set.
        5 Parameter 'a' and 't' should be set together.
HELP
}
sub VERSION_MESSAGE {
    print <<VERSION;
NEBPKM - Calculate BPKM of id-LncRNA seed regions' neighbor exons.
Version: $VERSION, Maintainer: $AUTHOR
VERSION
}
sub help {
    print <<SIMPLE_HELP;
Usage: NEBPKM.pl -r <reference file> -b <BAM file> -f <id-LncRNA file> -s <read number>
                [-e] [-c] [-i] [-a additional BAM file] [-t additional BAM file read number]
                [-l <read length>] [-d <output directory>] [-o <output file>]
See 'NEBPKM.pl --help' for more details.
SIMPLE_HELP
    exit;
}

# get options and check
my %opts = ( l=>'75' );
getopts('hecir:b:f:s:a:t:l:d:o', \%opts) or help();
help() if $opts{h} or not exist_essential_parameter();
if (exists $opts{a}) {
    croak "Parameter 'a' and 't' should be set together!" if not exists $opts{t};
}
if (exists $opts{t}) {
    croak "Parameter 'a' and 't' should be set together!" if not exists $opts{a};
}

# check and set output directory
if(exists $opts{d}) {
    $opts{d} =~ s/([^\/]$)/$1\//;
    -d $opts{d} or croak "The diresctor you input is wrong!";
}
else {
    $opts{d} = '';
}

# check and set output filename
my $basename = fileparse($opts{f});
$basename =~ /(.*).id-LncRNA/;
$opts{o} = $1 if not exists $opts{o};

# set output, log and error file
open my $result,'>',"$opts{d}$opts{o}.id-LncRNA_NE" or croak "Can't open file $opts{d}$opts{o}.id-LncRNA_NE!";
open STDOUT,'>',"$opts{d}$opts{o}_NE.log";
open STDERR,'>',"$opts{d}$opts{o}_NE.error";

# record command and starting time
my $command;
while (my ($key, $value) = each %opts) {
    $command .= " -$key $value";
}
print "Command: $0$command\n";
print 'Program starts at ', scalar(localtime(time)), "\n";

# read reference file and construct index
open my $ref_file,'<',$opts{r} or croak "Can't open $opts{r} file!";
my $refname = fileparse($opts{r}, qr/\Q.txt\E/);
my $gene_exon_sta_end_ref
    = ExtractInfo($ref_file, $refname, 'all', '-', 1, 1, $opts{e});
close $ref_file;

# record header
my @comp_header = ('', '', '', '', '', '');
if ($opts{c}) {
    @comp_header[1..5] = ("\tFoldChangeWithAdditionalBPKM","\tFoldChangeWithLeftExonBPKM", "\tFoldChangeWithAdditionalLeftExonBPKM", "\tFoldChangeWithRightExonBPKM", "\tFoldChangeWithAdditionalRightExonBPKM");
}
$comp_header[0] = "\tInteresting" if $opts{i};
$comp_header[0] .= "\tAdditionalInteresting" if $opts{a};
if (exists $opts{a}) {
    print $result "Position$comp_header[0]\tLength\tBPKM\tAdditionalBPKM$comp_header[1]\tLeftExon\tLeftExonInfo\tLeftExonLength\tLeftExonDistance\tLeftExonBPKM$comp_header[2]\tAdditionalLeftExonBPKM$comp_header[3]\tRightExon\tRightExonInfo\tRightExonLength\tRightExonDistance\tRightExonBPKM$comp_header[4]\tAdditionalRightExonBPKM$comp_header[5]\n";
}
else {
    print $result "Position$comp_header[0]\tLength\tBPKM\tLeftExon\tLeftExonInfo\tLeftExonLength\tLeftExonDistance\tLeftExonBPKM$comp_header[2]\tRightExon\tRightExonInfo\tRightExonLength\tRightExonDistance\tRightExonBPKM$comp_header[4]\n";
}

open my $ilrna_file,'<',$opts{f} or croak "Can't open $opts{f} file!";
while (<$ilrna_file>) { # for each record in id-LncRNA file

    next if $. == 1; # skip the header

    chomp;
    my @line = split /\s/,$_,4;
    my @gene = split /\s/,$line[3]; # get gene info
    my ($chr, $region) = split /:/,$line[0];

    my $exon_regions_ref = getexonregion($chr, \@gene, $gene_exon_sta_end_ref);
    my ($left_exon, $left_exon_info, $left_exon_length, $left_distance)
        = neighborexon($chr, $region, $exon_regions_ref, 0); # left exon
    my ($right_exon, $right_exon_info, $right_exon_length, $right_distance)
        = neighborexon($chr, $region, $exon_regions_ref, 1); # right exon

    my @foldchange = ('') x 5;
    my $left_bpkm = calculateBPKM($left_exon, $left_exon_length, $opts{b}, $opts{s}); # left BPKM
    my $right_bpkm = calculateBPKM($right_exon, $right_exon_length, $opts{b}, $opts{s}); # right BPKM
    if ($opts{c}) { # calculate fold change
        $foldchange[0] = compare($line[2], $left_bpkm);
        $foldchange[1] = compare($line[2], $right_bpkm);
    }
    my ($additional_bpkm, $additional_left_bpkm, $additional_right_bpkm);
    if (exists $opts{a}) {
        $additional_bpkm = calculateBPKM($line[0], $line[1], $opts{a}, $opts{t}); # additional BPKM
        $additional_left_bpkm = calculateBPKM($left_exon, $left_exon_length, $opts{a}, $opts{t}); # left additional BPKM
        $additional_right_bpkm = calculateBPKM($right_exon, $right_exon_length, $opts{a}, $opts{t}); # right additional BPKM
        if ($opts{c}) { # calculate fold change
            $foldchange[2] = compare($line[2], $additional_bpkm);
            $foldchange[3] = compare($line[2], $additional_left_bpkm);
            $foldchange[4] = compare($line[2], $additional_right_bpkm);
        }
    }

    # all interesting flag
    my $interesting_flag = '';
    if ($opts{i}) {
        # interesting flags
        my @interesting = map {($_ eq '###' or $_ >= 2) ? 1 : 0} @foldchange[0, 1];
        $interesting_flag = sum(@interesting) == 2 ? "\t1" : "\t0" if $opts{c};
        if ($interesting_flag eq "\t1") {
            $interesting_flag .= ($foldchange[2] eq '###' or $foldchange[2] >= 2) ? "\t1" : "\t0" if $opts{a};
        }
        else {
            $interesting_flag .= "\t0" if $opts{a};
        }
    }

    @foldchange = map {"\t$_"} @foldchange if $opts{c};

    # record result
    if (exists $opts{a}) {
        print $result "$line[0]$interesting_flag\t$line[1]\t$line[2]\t$additional_bpkm$foldchange[2]\t$left_exon\t$left_exon_info\t$left_exon_length\t$left_distance\t$left_bpkm$foldchange[0]\t$additional_left_bpkm$foldchange[3]\t$right_exon\t$right_exon_info\t$right_exon_length\t$right_distance\t$right_bpkm$foldchange[1]\t$additional_right_bpkm$foldchange[4]\n";
    }
    else {
        print $result "$line[0]$interesting_flag\t$line[1]\t$line[2]\t$left_exon\t$left_exon_info\t$left_exon_length\t$left_distance\t$left_bpkm$foldchange[0]\t$right_exon\t$right_exon_info\t$right_exon_length\t$right_distance\t$right_bpkm$foldchange[1]\n";
    }

}

# close file handle
close $ilrna_file;
close $result;

# delete error file if no errors
unlink "$opts{d}$opts{o}_NE.error" if -z "$opts{d}$opts{o}_NE.error";

# record ending time
print 'Program ends at ', scalar(localtime(time)), "\n";

##########Subroutine##########

#
# Function: Check essential parameters.
#
sub exist_essential_parameter {
    my @essential_parameter = qw( r b f s );
    exists $opts{$_} || return 0 for @essential_parameter;
    return 1;
}

#
# Function: Get exon region
#
sub getexonregion {
    my ($chr, $gene_ref, $gene_exon_edge_ref) = @_;
    my @gene_index;
    push @gene_index,$chr . ':' . $_ for @{$gene_ref};
    my @exon_regions;
    for (@gene_index) {
        s/-/#/g;
        if (exists $gene_exon_edge_ref->{$_}) {
            push @exon_regions,@{$gene_exon_edge_ref->{$_}};
        }
    }
    return \@exon_regions;
}

#
# Function: Extract neighbor exon region
#
sub neighborexon {
    my ($chr, $region, $exon_regions_ref, $flag) = @_;
    my @sorted_exon_regions;
    if ($flag) {
        @sorted_exon_regions = map { $_->[0] }
                               sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
                               map { [ $_, split /-/ ] }
                               @$exon_regions_ref;
    }
    else {
        @sorted_exon_regions = map { $_->[0] }
                               sort { $b->[2] <=> $a->[2] or $b->[1] <=> $a->[1] }
                               map { [ $_, split /-/ ] }
                               @$exon_regions_ref;
    }
    my ($left, $right) = split /-/,$region;
    my ($exon, $gene_info, $length, $distance) = ('None') x 4;
    for (@sorted_exon_regions) {
        my ($left_edge, $right_edge, $info) = split /-/;
        $info =~ s/#/-/g;
        if ($flag) {
            if ($right <= $left_edge) {
                $exon = $chr . ':' . $left_edge . '-' . $right_edge;
                $gene_info = $info;
                $length = $right_edge - $left_edge;
                $distance = $left_edge - $right;
                last;
            }
        }
        else {
            if ($left >= $right_edge) {
                $exon = $chr . ':' . $left_edge . '-' . $right_edge;
                $gene_info = $info;
                $length = $right_edge - $left_edge;
                $distance = $left - $right_edge;
                last;
            }
        }
    }
    return ($exon, $gene_info, $length, $distance);
}

#
# Function: Calculate BPKM
#
sub calculateBPKM {
    my $region = shift;
    my $length = shift;
    my $bam_file = shift;
    my $size = shift;
    return 0 if $region eq 'None';
    open my $read_iterator,'-|',"samtools view $bam_file $region";
    my @read_segments;
    while(<$read_iterator>) {
        my ($pos, $cigar) = (split)[3, 5];
        my $segments_ref = ReadSplit($pos, $cigar, '-');
        push @read_segments,@{$segments_ref};
    }
    close $read_iterator;
    return 0 if not exists $read_segments[0];
    $region =~ s/^chr(\d+|X|Y|M)://;
    my $reads_in_region_ref = OverlapMap([$region], \@read_segments, '-');
    my $base = 0;
    for (@{$reads_in_region_ref}) {
        my ($sta, $end) = split /-/;
        $base += $end - $sta;
    }
    return ($base * 10**9) / ($size * $opts{l} * $length);
}

#
# Function: Compare BPKMs to find interesting seeds
#
sub compare {
    my ($bpkm1, $bpkm2) = @_;
    if ($bpkm2 == 0) {
        return '###';
    }
    else {
        return ($bpkm1 / $bpkm2);
    }
}
