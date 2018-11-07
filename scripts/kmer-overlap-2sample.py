#!/usr/bin/env python
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,no-member
from __future__ import print_function
from __future__ import absolute_import

import argparse
import sys
import os

import screed
from screed.fasta import fasta_iter

import khmer
from khmer import khmer_args
from khmer.khmer_args import (build_nodegraph_args, DEFAULT_MAX_TABLESIZE)
from khmer.utils import broken_paired_reader, write_record

def kmer_degree(kmer, ht):
    s = "ACGT"
    left = sum(ht.get('{}{}'.format(i, kmer)) for i in s)
    right = sum(ht.get('{}{}'.format(kmer, i)) for i in s)
    return left, right

def main():
    parser = build_nodegraph_args()
    parser.add_argument('readfile1', help='fasta sequence file to be loaded in hashtable, use "-" if from stdin')
    parser.add_argument('readfile2', help='fasta readfile to query against hashtable, use "-" if from stdin')
    parser.add_argument('--shared', help='shared kmer in readfile 1 and 2')
    parser.add_argument('--uniq2', help='uniq kmer in readfile2')
    parser.add_argument('--x2', default='1e8', help='max_table size for readfile2')
    parser.add_argument('--N2', default='4', help='# of table (N) for readfile2')

    args = parser.parse_args()
    print(args)

    K = args.ksize
    HT_SIZE = args.max_tablesize
    N_HT = args.n_tables
    HT_SIZE2 = int(float(args.x2))
    N_HT2 = int(args.N2)

    # positional
    readfile1 = args.readfile1
    readfile2 = args.readfile2
    shared = args.shared
    uniq2 = args.uniq2
    if readfile1 == '-' and readfile2 == '-':
        mes = ('*** Only one of readfile1 and readfile2 '
                 'can be read from stdin')
        print(mes, file=sys.stderr)

    try:
        if readfile1 == '-':
            fp1 = sys.stdin
        else:
            fp1 = open(readfile1)

        if readfile2 == '-':
            fp2 = sys.stdin
        else:
            fp2 = open(readfile2)

        if uniq2:
            fw2 = open(uniq2, 'w')

        # create a hashbits data structure
        ht = khmer.Nodetable(K, HT_SIZE, N_HT)

        # load contigs, connect into N partitions
        print('loading input reads from {}..'.format(os.path.basename(readfile1)),
                                    file=sys.stderr)
        #ht.consume_seqfile(readfile1)
        for record in fasta_iter(fp1):
            ht.consume(record.sequence)

        # Change 0.2 only if you really grok it.  HINT: You don't.
        fp_rate = khmer.calc_expected_collisions(ht)
        print('fp rate estimated to be {:1.3f}'.format(fp_rate), file=sys.stderr)

        if fp_rate > 0.01:
            mes = ('**\n'
                   '** ERROR: the counting hash is too small for\n'
                   '** {}.  Increase hashsize/num ht.\n'
                   '**\n'
                   '** Do not use these results!!')
            print(mes.format(os.path.basename(readfile1)), file=sys.stderr)
            sys.exit(-1)


        n_unique1 = ht.n_unique_kmers()
        # create a hashbits data structure
        ht2 = khmer.Nodetable(K, HT_SIZE2, N_HT2)

        n_unique2 = 0
        n_shared = 0

        for n, record in enumerate(fasta_iter(fp2)):
            name = record['name']
            sequence = record['sequence']
            seq_len = len(sequence)
            for i in range(0, seq_len + 1 - K):
                kmer = sequence[i:i + K]

                if (not ht2.get(kmer)):
                    n_unique2 += 1
                    if ht.get(kmer):
                        n_shared += 1
                    else:
                        mes = '>{}__{}  length_{};k_{}\n{}\n'
                        fw2.write(mes.format(name, i, seq_len, K, kmer))
                ht2.count(kmer)

        mes = ('Unique kmer in {}:\t{}\n'
               'Shared kmer:\t{}\n'
               'Unique kmer in {}:\t{}\n')

        print(mes.format(os.path.basename(readfile1), n_unique1, 
                           n_shared, 
                         os.path.basename(readfile2), n_unique2))

    finally:
        fp1.close()
        fp2.close()
        fw2.close()

if __name__ == '__main__':
    main()
