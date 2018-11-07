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
    parser.add_argument('query', help=('fasta readfile to query against' 
                                         'hashtable, use "-" if from stdin'))
    parser.add_argument('ref', nargs='+', 
                          help='fasta sequence file to be loaded in hashtable')
    parser.add_argument('--x2', default='1e8', 
                          help='max_table size for readfile2')
    parser.add_argument('--N2', default='4', 
			  help='# of table (N) for readfile2')

    args = parser.parse_args()
    #print(args, file=sys.stderr)

    K = args.ksize
    HT_SIZE = args.max_tablesize
    N_HT = args.n_tables
    HT_SIZE2 = int(float(args.x2))
    N_HT2 = int(args.N2)

    # positional
    query = args.query
    refs = args.ref
    print('{} refs to be loaded'.format(len(refs)), file=sys.stderr)
    if query == '-' and refs == ['-']:
        print('*** query and ref can not both be "-" (read from stdin)', 
                file=sys.stderr)
    # create a hashbits data structure
    start_time = time.time()
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

    pair = 0
    forward = 0
    reverse = 0
    other = 0
    total_pair = 0
    for n, is_pair, r1, r2 in broken_paired_reader(khmer.ReadParser(query, 
                                                        require_paired=True)):
    #for n, record in enumerate(screed.open(query)):
	total_pair += 1
        share_list = []
        for record in [r1, r2]:
            name, desc = record.name.split(None, 1)
            sequence = record.sequence.replace('N', 'A')
            seq_len = len(sequence)
            if seq_len < K:
                print('*** {} is shorter than {}..'.format(r1.name, K), 
                        file=sys.stderr)
                continue
            for i in range(0, seq_len + 1 - K):
                kmer = sequence[i:i + K]
                if ht.get(kmer):
                    share_list.append(1)
                    break
                else:
                    share_list.append(0)

        if share_list == [1, 1]:
            pair += 1 
        elif share_list == [1, 0]:
            forward += 1 
        elif share_list == [0, 1]:
            reverse += 1 
        else: #[0, 0]
            other += 1
            # do not print
            continue

        mes = ('>{}  {}||uniq_{}\n{}\n'
                '>{}  {}||uniq_{}\n{}')
        l1 = r1.name.split(None, 1)
        l2 = r2.name.split(None, 1)
        print(mes.format(l1[0], l1[1], share_list[0], r1.sequence,
                            l2[0], l2[1], share_list[1], r2.sequence))


    mes = ('Unique kmer in ref:\t{}\n'
           'Total pair:\t{}\n'
           'Both primers uniq:\t{}\n'
           'Pair with forward uniq:\t{}\n'
           'Pair with reverse uniq:\t{}')

    print(mes.format(n_unique1, total_pair, pair, forward, reverse), 
                        file=sys.stderr)


if __name__ == '__main__':
    main()
