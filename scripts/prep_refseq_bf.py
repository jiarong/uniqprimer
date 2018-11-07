#/usr/bin/env python
# tested on khmer/2.1.1, screed/1.0
from __future__ import print_function, unicode_literals
from __future__ import absolute_import

import argparse
import sys
import os
import json
import threading
import textwrap
import time

import screed
from screed.fasta import fasta_iter

import khmer
from khmer import khmer_args
from khmer.khmer_args import (build_nodegraph_args, report_on_config, info,
                              add_threading_args, calculate_graphsize,
                              sanitize_help, DEFAULT_MAX_TABLESIZE)
from khmer.utils import broken_paired_reader, write_record, clean_input_reads
from khmer.kfile import check_file_writable
from khmer.kfile import check_input_files
from khmer.kfile import check_space_for_graph

def kmer_degree(kmer, ht):
    s = "ACGT"
    left = sum(ht.get('{}{}'.format(i, kmer)) for i in s)
    right = sum(ht.get('{}{}'.format(kmer, i)) for i in s)
    return left, right

def main():
    parser = build_nodegraph_args("find uniq kmer in query compard to refs")
    parser.add_argument('ref', nargs='+', 
                     help='fasta sequence file to be loaded in bloom filter')
    parser.add_argument('--bfout', default='nodetable.bf', 
                     help='output bloom filter of ref')

    args = parser.parse_args()

    K = args.ksize
    HT_SIZE = args.max_tablesize
    N_HT = args.n_tables

    # positional
    refs = args.ref

    start_time = time.time()
    print('{} refs to be loaded'.format(len(refs)), file=sys.stderr)
    ht = khmer.Nodetable(K, HT_SIZE, N_HT)
    end_time = time.time()
    secs = end_time - start_time
    mes = 'initiation of bloom filter took {:.2f} hours..'
    print(mes.format(secs/3600.0), file=sys.stderr)

    for index, filename in enumerate(refs):
        if index != 0 and index % 100 == 0:
            end_time = time.time()
            secs = end_time - start_time
            mes = '{} refs have been loaded with in {:.2f} hours ..'
            print(mes.format(index, secs/3600.0), file=sys.stderr)
        try:
            ht.consume_seqfile(filename)
        except OSError as e:
            mes = ('*** Skipping due to OSError (machine or system problem):'
                       ' {}\n'
                   '*** Detailed error message:\n'
                   '*** {}')
            print(mes.format(os.path.basename(filename), str(e)), 
                                                    file=sys.stderr)
            continue

    print('Saving bloom filter to {}..'.format(args.bfout), 
            file=sys.stderr)
    ht.save(args.bfout)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    mes = 'fp rate estimated to be {:1.3f}'
    print(mes.format(fp_rate), file=sys.stderr)

    if fp_rate > 0.01:
        mes = ('**\n'
               '** ERROR: the counting hash is too small for\n'
               '** refs.  Increase hashsize/num ht.\n'
               '**\n'
               '** Do not use these results!!')
        sys.exit(-1)

    n_unique1 = ht.n_unique_kmers()

    mes = ('Unique kmer:\t{}\n')

    print(mes.format(n_unique1), file=sys.stderr)


if __name__ == '__main__':
    main()
