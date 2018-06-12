#!/usr/bin/env python3

'''
title: SizeModify.py
description: Modify peak size.
author: Xiao-Ou Zhang
version: 0.2.0
'''

import sys
sys.path.insert(0, '/picb/rnomics1/xiaoou/program/usefullib/python/')
sys.path.insert(0, '/picb/rnomics1/xiaoou/program/usefullib/test/')
from interval import Interval
from map import overlapwith
from bpkm import calculatebpkm
import argparse
import os
import os.path
import re


def argumentparse():
    '''
    argumentparse() -> parser -- create a parser
    '''
    parser = argparse.ArgumentParser(description='Modify peak size.')
    parser.add_argument('-v', action='version', version='%(prog)s 0.2.0')
    parser.add_argument('-p', required=True, metavar='Peak File',
                        dest='peakfile')
    parser.add_argument('-r', required=True, metavar='Reference File',
                        dest='referencefile')
    parser.add_argument('-b', required=True, metavar='BAM File',
                        dest='bamfile')
    parser.add_argument('-t', required=True, metavar='Read Number',
                        dest='size')
    parser.add_argument('-m', required=True,
                        metavar='Multiple Mapping BAM File', dest='mbamfile')
    parser.add_argument('-l', metavar='Read Length', default=75, type=int,
                        help='(default 75)', dest='length')
    parser.add_argument('-n', metavar='Chromosome Number', default='24',
                        help='(default 24)', dest='number')
    parser.add_argument('-d', metavar='Out File Directory',
                        dest='outfiledirectory')
    return parser


def extractminintron(rfile, file_type, num, xy):
    '''
    extractminintron(rfile, file_type, num, xy) -> minintron
    -- extract min intron region
    '''
    minintron = {'chr' + str(i): [] for i in list(range(1, num)) + xy}
    if file_type == 'knownGene.txt':
        exon, isoform = refparse(rfile, 0, num, xy)
    elif file_type == 'refFlat.txt':
        exon, isoform = refparse(rfile, 1, num, xy)
    for i in list(range(1, num)) + xy:
        exon_interval = Interval(exon['chr{}'.format(i)])
        isoform_interval = Interval(isoform['chr{}'.format(i)])
        intron_interval = isoform_interval - exon_interval
        for interval in intron_interval.interval:
            interval.append('chr{0}:{1[0]}-{1[1]}'.format(i, interval))
        minintron['chr{}'.format(i)] = intron_interval.interval
    return minintron


def refparse(rfile, flag, num, xy):
    '''
    refparse(rfile, flag, num, xy) -> (exon, isoform) -- parse reference file
    '''
    exon = {'chr{}'.format(i): [] for i in list(range(1, num)) + xy}
    isoform = {'chr{}'.format(i): [] for i in list(range(1, num)) + xy}
    chromosome = ['chr{}'.format(i) for i in list(range(1, num)) + xy]
    with open(rfile, 'r') as ref:
        for line in ref:
            chrom = line.split()[1 + flag]
            if chrom not in chromosome:
                continue
            sta, end = line.split()[(3 + flag):(5 + flag)]
            isoform[chrom].append([sta, end])
            exonsta = line.split()[8 + flag].split(',')[:-1]
            exonend = line.split()[9 + flag].split(',')[:-1]
            for s, e in zip(exonsta, exonend):
                exon[chrom].append([s, e])
    return (exon, isoform)


def statistics(data, chrom, out):
    '''
    statistics(data, chrom, out) -> result -- statistics
    '''
    output_buffer = ''
    result = []
    out_buffer = []
    tmp = []
    tmp_id = set([])
    n = 0
    write_flag = 0
    while True:
        try:
            p = data[n]
        except IndexError:
            if write_flag:
                if out_buffer:
                    result.append(out_buffer + [tmp[3]])
            else:
                out.write(output_buffer)
            break
        if len(p) <= 4:
            intron_id = set([])
            intron_region = ''
        else:
            intron_id = set(p[4:])
            intron_region = '({})'.format(p[-1])
        if intron_id & tmp_id == set([]):
            if write_flag:
                result.append(out_buffer + [tmp[3]])
                write_flag = 0
            else:
                out.write(output_buffer)
            out_buffer = [intron_region]
            output_buffer = ''
        else:
            distance = p[0] - tmp[1]
            out_buffer.append(str(distance))
            write_flag = 1
        out_buffer.extend(['{}:{}-{}'.format(chrom, p[0], p[1]), p[2]])
        length = p[1] - p[0]
        output_buffer += '{0}:{1[0]}-{1[1]}\t{2}\t{1[2]}\t{1[3]}\n'.format(
                chrom, p, length)
        tmp = p
        tmp_id = intron_id
        n += 1
    return result


def combine(data, chrom, bam, mbam, size, length):
    '''
    combine(data, chrom, bam, mbam, size, length)
    -> (combined_peak, cpeak, rpeak) -- combine peak
    '''
    combined_peak = []
    cpeak = []
    rpeak = []
    for region in data:
        record = []
        for i in range(1, int((len(region) - 4) / 3) + 1):
            left = region[(3 * i - 2):(3 * i)]
            right = region[(3 * i + 1):(3 * i + 3)]
            sta, left_edge = re.split(r':|-', region[3 * i - 2])[1:3]
            right_edge, end = re.split(r':|-', region[3 * i + 1])[1:3]
            tmp = calculatebpkm(chrom, left_edge, right_edge, mbam,
                                        getsegment=True)
            if tmp != 0:
                tmp = Interval(tmp)
            else:
                tmp = Interval([])
            if ((tmp.interval == [[int(left_edge), int(right_edge)]] or
                    left_edge == right_edge) and
                    (1 / 3 <= float(left[1]) / float(right[1]) <= 3)):
                record.append([sta, end])
                cpeak.append([region[0]] + left + [region[3 * i]] + right)
            else:
                record.extend([[sta, left_edge], [right_edge, end]])
                rpeak.append([region[0]] + left + [region[3 * i]] + right)
        record = Interval(record)
        for interval in record.interval:
            bpkm = calculatebpkm(chrom, interval[0], interval[1], bam, size,
                                length)
            combined_peak.append(['{0}:{1[0]}-{1[1]}'.format(chrom, interval),
                                    str(interval[1] - interval[0]), str(bpkm),
                                    region[-1]])
    return (combined_peak, cpeak, rpeak)


def output(data, f):
    '''
    output(data, f) -> None -- output result
    '''
    for item in data:
        f.write('\t'.join(item) + '\n')


if __name__ == '__main__':
    parser = argumentparse()
    opts = vars(parser.parse_args(sys.argv[1:]))
    # set output directory and file name
    if opts['outfiledirectory'] is None:
        opts['outfiledirectory'] = os.getcwd()
    else:
        if not os.path.isdir(opts['outfiledirectory']):
            print('Your out file directory is wrong!')
            sys.exit(1)
    opts['outfiledirectory'] = re.sub(r'/?$', r'/', opts['outfiledirectory'])
    file_type = os.path.split(opts['referencefile'])[1]
    if file_type not in ('refFlat.txt',
            'knownGene.txt') or not os.path.isfile(opts['referencefile']):
        print('The reference file is wrong!')
        sys.exit(1)
    name, suffix = os.path.splitext(os.path.split(opts['peakfile'])[1])
    outf1 = open(opts['outfiledirectory'] + name + '_size_statistics_use.txt',
                'w')
    outf2 = open(opts['outfiledirectory'] + name + '_size_statistics_rest.txt',
                'w')
    outf = open(opts['outfiledirectory'] + name + '_sizemodified' + suffix,
                'w')
    if opts['number'].endswith('X'):
        xy = ['Y']
        num = int(opts['number'][:-1])
    elif opts['number'].endswith('Y'):
        xy = ['X']
        num = int(opts['number'][:-1])
    else:
        xy = ['X', 'Y']
        num = int(opts['number']) - 1
    # get min intron region
    minintron = extractminintron(opts['referencefile'], file_type, num, xy)
    # parse peak file
    peak = {'chr{}'.format(i): [] for i in list(range(1, num)) + xy}
    if not os.path.isfile(opts['peakfile']):
        print('Please input the right peak file!')
        sys.exit(1)
    else:
        with open(opts['peakfile'], 'r') as pfile:
            outf.write(pfile.readline())
            for line in pfile:
                chrom, sta, end = re.split(r':|-', line.split()[0])
                bpkm, info = line.split('\t')[2:4]
                info = info.rstrip()
                peak[chrom].append([sta, end, bpkm, info])
    # annotate peak
    annotated_peak = {'chr{}' + str(i): [] for i in list(range(1, num)) + xy}
    classified_peak = {'chr{}' + str(i): [] for i in list(range(1, num)) + xy}
    for i in list(range(1, num)) + xy:
        annotated_peak['chr' + str(i)] = overlapwith(peak['chr' + str(i)],
                                                    minintron['chr' + str(i)])
        chrom = 'chr' + str(i)
        classified_peak[chrom] = statistics(annotated_peak[chrom], chrom, outf)
        combined_peak, cpeak, rpeak = combine(classified_peak[chrom], chrom,
                                            opts['bamfile'], opts['mbamfile'],
                                            opts['size'], opts['length'])
        output(cpeak, outf1)
        output(rpeak, outf2)
        output(combined_peak, outf)
    outf1.close()
    outf2.close()
    outf.close()
