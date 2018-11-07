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

    parser.add_argument('--x2', default='1e8', 
                          help='max_table size for readfile2')
    parser.add_argument('--N2', default='4', 
                            help='# of table (N) for readfile2')

    parser.add_argument('--bfout', help='output bloom filter of ref')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--shared', dest='output', action='store_const', 
                        const='shared', help='output shared kmers')
    group.add_argument('--uniq', dest='output', action='store_const', 
                        const='uniq', help='output uniq kmers in query')

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument('--ref', nargs='+', 
                        help='fasta sequence file to be loaded in bloom filter')
    group2.add_argument('--load', help='load existing bloom filter')

    parser.set_defaults(output='uniq')
    args = parser.parse_args()
    #print(args, file=sys.stderr)

    K = args.ksize
    HT_SIZE = args.max_tablesize
    N_HT = args.n_tables
    HT_SIZE2 = int(float(args.x2))
    N_HT2 = int(args.N2)

    # positional
    query = args.query
    output = args.output

    start_time = time.time()
    # load from existing bloom filter
    if args.load:
        print('loading bloom filter from {}..'.format(args.load),
                file=sys.stderr)
        ht = khmer.load_nodetable(args.load)
        k = ht.ksize()
        mes = ('*** incompatible ksize ({}) in {} with parameters K on '
                'command line ({})')
        assert k == K, mes.format(k, args.load, K)
        end_time = time.time()
        secs = end_time - start_time
        mes = 'load bloom filter ({}) took {:.2f} hours..'
        print(mes.format(os.path.basename(args.load), secs/3600.0), 
                file=sys.stderr)

    # create a hashbits data structure
    else:
        refs = args.ref
        print('{} refs to be loaded'.format(len(refs)), file=sys.stderr)
        if query == '-' and refs == ['-']:
            print('*** query and ref can not both be "-" (read from stdin)', 
                    file=sys.stderr)
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

        if args.bfout:
            if args.load:
                mes = '*** Bloom filter exists as {}, NOT saving again as {}..'
                print(mes.format(args.load, args.bfout), file=sys.stderr)
            else:
                print('*** Saving bloom filter to {}..'.format(args.bfout), 
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


    # create a hashbits data structure
    ht2 = khmer.Nodetable(K, HT_SIZE2, N_HT2)

    n_unique2 = 0
    n_shared = 0

    if output == 'uniq':
        for n, record in enumerate(khmer.ReadParser(query)):
        #for n, record in enumerate(screed.open(query)):
            _l = record.name.split(None, 1)
            if len(_l) == 2:
                name, desc = _l
            else:
                name = _l[0]
                desc = ''
            sequence = record.sequence.replace('N', 'A')
            seq_len = len(sequence)
            if seq_len < K:
                continue
            for i in range(0, seq_len + 1 - K):
                kmer = sequence[i:i + K]

                if (not ht2.get(kmer)):
                    n_unique2 += 1
                    if ht.get(kmer):
                        n_shared += 1
                    else:
                        mes = '>{}__{}  {}||length_{};k_{}\n{}'
                        print(mes.format(name, i, desc, seq_len, K, kmer))
                ht2.count(kmer)

    elif output == 'shared':
        for n, record in enumerate(khmer.ReadParser(query)):
        #for n, record in enumerate(screed.open(query)):
            _l = record.name.split(None, 1)
            if len(_l) == 2:
                name, desc = _l
            else:
                name = _l[0]
                desc = ''
            sequence = record.sequence.replace('N', 'A')
            seq_len = len(sequence)
            if seq_len < K:
                continue
            for i in range(0, seq_len + 1 - K):
                kmer = sequence[i:i + K]

                if (not ht2.get(kmer)):
                    n_unique2 += 1
                    if ht.get(kmer):
                        n_shared += 1
                        mes = '>{}__{}  {}||length_{};k_{}\n{}'
                        print(mes.format(name, i, desc, seq_len, K, kmer))
                    else:
                        pass

                ht2.count(kmer)

    mes = ('Unique kmer in {} (query):\t{}\n'
           'Shared kmer:\t{}\n'
           'Unique kmer in {}:\t{}\n')

    print(mes.format(os.path.basename(query), n_unique2, 
                       n_shared, 
                     'refs', n_unique1),
                       file=sys.stderr)


if __name__ == '__main__':
    main()
