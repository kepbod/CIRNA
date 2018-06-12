#!/usr/bin/env python3

'''
title: AnnotateRNA.py
description: Annotate RNA info for each seed.
author: Xiao-Ou Zhang
version: 1.0.0
'''

import argparse
import sys
import os.path
import os
import re
sys.path.insert(0, '/picb/rnomics1/xiaoou/program/usefullib/python/')
from map import overlapwith


def argumentparse():
    '''
    argumentparse() -> parser -- create a parser
    '''
    parser = argparse.ArgumentParser(
            description='Annotate RNA info for each seed.')
    parser.add_argument('-v', action='version', version='%(prog)s 1.0.0')
    parser.add_argument('-r', required=True, metavar='RNA Info File',
                        dest='rnainfo')
    parser.add_argument('-s', required=True, metavar='Seed File',
                        dest='seedfile')
    parser.add_argument('-e', action='store_true', dest='extending',
                        help='extending flag (default false)')
    parser.add_argument('-n', metavar='Chromosome Size', default=24, type=int,
                        help='(default 24)', dest='chromosomesize')
    parser.add_argument('-l', metavar='Extending Size', default=500, type=int,
                        help='(default 500)', dest='extendingsize')
    parser.add_argument('-d', metavar='Out File Directory',
                        dest='outfiledirectory')
    parser.add_argument('-o', metavar='Out File Name', dest='outfilename')
    return parser


def getrnainfofile(rnainfofile):
    '''
    getrnainfofile(rnainfofile) -> filehandle -- check the rna info file and
    return file handle
    '''
    if not os.path.isfile(rnainfofile):
        print('Please input the right rna info file!')
        sys.exit(1)
    else:
        return open(rnainfofile, 'r')


def parsernainfofile(rnainfo, file_type, size):
    '''
    parsernainfofile(rnainfo, file_type, size) -> rna_list -- parse rna info
    file and extract RNA list
    '''
    rna_list = {'chr' + str(i): [] for i in range(1, size - 1)}
    rna_list['chrX'] = []
    rna_list['chrY'] = []
    # extract RNA info according to chromosome info
    for rna in rnainfo:
        if file_type == 'wgRna.txt':
            chrom, sta, end, name = rna.split()[1:5]
            rnatype = rna.split()[-1]
        elif file_type == 'rmsk.txt':
            chrom, sta, end = rna.split()[5:8]
            name = rna.split()[10]
            rnatype = ':'.join(rna.split()[11:13])
        else:
            chrom, sta, end, name = rna.split()[0:4]
            rnatype = 'TF'
        if chrom not in rna_list:
            continue
        rna_list[chrom].append([sta, end, ':'.join([name, rnatype])])
    return rna_list


def getseedfile(seed_file):
    '''
    getseedfile(seed_file) -> filehandle -- check the seed file and return file
    handle
    '''
    if not os.path.isfile(seed_file):
        print('Please input the right seed file!')
        sys.exit(1)
    else:
        return open(seed_file, 'r')


def parseseedfile(seed_file, size, extending_flag, extending_size):
    '''
    parseseedfile(seed_file, size, extending_flag, extending_size)
    ->  seed_list, extending_seed_list, header  -- parse seed file and extract
    seed list
    '''
    header = seed_file.readline().rstrip()  # extract head in case write back
    seed_list = {'chr' + str(i): [] for i in range(1, size - 1)}
    seed_list['chrX'] = []
    seed_list['chrY'] = []
    extending_seed_list = {'chr' + str(i): [] for i in range(1, size - 1)}
    extending_seed_list['chrX'] = []
    extending_seed_list['chrY'] = []
    for seed in seed_file:  # check each seed
        seed = seed.rstrip()
        chrom, sta_end = seed.split()[0].split(":")
        sta, end = sta_end.split('-')
        extending_sta = int(sta) - extending_size
        extending_end = int(end) + extending_size
        seed_list[chrom].append([sta, end, seed])
        if extending_flag:
            extending_seed_list[chrom].append([extending_sta, extending_end])
    return seed_list, extending_seed_list, header


def overlap(seed_list, rna_list, size):
    '''
    overlap(seed_list, rna_list, size) -> None
    -- check overlaps between seeds and RNA regions
    '''
    for chrom in list(range(1, size - 1)) + ['X', 'Y']:
        index = seed_list['chr' + str(chrom)]
        interval = rna_list['chr' + str(chrom)]
        seed_list['chr' + str(chrom)] = overlapwith(index, interval)


def writeoutput(lst, elt, header, out_directory, out_name, file_type,
                extending_flag, size):
    '''
    writeoutput(output_list, extending_output_list, header, out_directory,
    out_name, file_type, extending_flag, size) -> None -- write outputs
    '''
    out_directory = re.sub(r'/?$', r'/', out_directory)
    with open(out_directory + out_name, 'w') as outf:
        if file_type == 'wgRna.txt':
            added_header = '\tRNAInfo'
            extend_header = '\tExtendRNAInfo'
        elif file_type == 'rmsk.txt':
            added_header = '\tRepeats'
            extend_header = '\tExtendRepeatInfo'
        else:
            added_header = '\tTFInfo'
            extend_header = '\tExtendTFInfo'
        added_header += extend_header if extending_flag else ''
        outf.write(header + added_header + '\n')
        for chrom in list(range(1, size - 1)) + ['X', 'Y']:
            index = 'chr' + str(chrom)
            if extending_flag:
                for a, b in zip(lst[index], elt[index]):
                    if a[3:]:
                        outf.write('{}\t|{}|'.format(a[2], '|'.join(a[3:])))
                    else:
                        outf.write('{}\t{}'.format(a[2], 'None'))
                    if b[2:]:
                        outf.write('\t|{}|'.format('|'.join(b[2:])))
                    outf.write('\n')
            else:
                for a in lst[index]:
                    if a[3:]:
                        outf.write('{}\t|{}|\n'.format(a[2], '|'.join(a[3:])))
                    else:
                        outf.write('{}\t{}\n'.format(a[2], 'None'))


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
    file_type = os.path.split(opts['rnainfo'])[1]
    if file_type not in ('wgRna.txt', 'rmsk.txt',
                        'wgEncodeRegTfbsClusteredV2.bed'):
        print('The rna info file is wrong!')
        sys.exit(1)
    if opts['outfilename'] is None:
        name, suffix = os.path.splitext(os.path.split(opts['seedfile'])[1])
        if file_type == 'wgRna.txt':
            opts['outfilename'] = name + '_rnainfo' + suffix
        elif file_type == 'rmsk.txt':
            opts['outfilename'] = name + '_repeatinfo' + suffix
        else:
            opts['outfilename'] = name + '_tfinfo' + suffix
    rnainfo = getrnainfofile(opts['rnainfo'])
    rna_list = parsernainfofile(rnainfo, file_type, opts['chromosomesize'])
    rnainfo.close()
    seed_file = getseedfile(opts['seedfile'])
    seed_list, extending_seed_list, header = parseseedfile(
            seed_file, opts['chromosomesize'], opts['extending'],
            opts['extendingsize'])
    seed_file.close()
    if opts['extending']:
        overlap(extending_seed_list, rna_list, opts['chromosomesize'])
    overlap(seed_list, rna_list, opts['chromosomesize'])
    writeoutput(seed_list, extending_seed_list, header,
                opts['outfiledirectory'], opts['outfilename'], file_type,
                opts['extending'], opts['chromosomesize'])
