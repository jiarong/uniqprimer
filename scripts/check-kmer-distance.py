#!/usr/bin/env python
from __future__ import print_function

import sys
import os
import numpy

import screed

def process(name, l, lower=150, upper=1000):
    l = sorted(set(l))
    arr = numpy.array(l)
    for i, n in enumerate(arr):
        distances = arr[i+1:] - n 
        for j, m in enumerate(distances):
            if m in range(lower, upper):
                #yield name, n, arr[i+1+j]
                print('{}\t{}\t{}'.format(name, n, arr[i+1+j]))

def main():
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} file.uniq2ref'
        print(mes.format(os.path.basename(sys.argv[0])), file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    lower = int(sys.argv[2])
    upper = int(sys.argv[3])
    for rec in screed.open(infile):
        name = rec.name.split(None, 1)[0]  #1095362_contig_14__2939
        _l = name.split('__')
        if len(_l) < 2:
            print('*** {} does have wrong format, skipping..'.format(name), 
                    file=sys.stderr)
            continue
        contig_name = _l[0]
        pos = int(_l[1])
        try:
            if contig_name == prev_contig_name:
                pos_list.append(pos)
            else:
                process(prev_contig_name, pos_list, lower, upper)
                # initialization 
                del pos_list
                pos_list = [pos]
                prev_contig_name = contig_name

        except UnboundLocalError as e:
            # for first record
            #UnboundLocalError due to prev_contig_name not defined
            pos_list = [pos]
            prev_contig_name = contig_name

    # process last batch
    try:
        process(prev_contig_name, pos_list, lower, upper)
    except UnboundLocalError as e:
        print('*** primer seqfile is empty, skipping..', file=sys.stderr)


if __name__ == '__main__':
    main()
