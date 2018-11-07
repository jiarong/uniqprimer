#! /usr/bin/env python
# works with screed >= 0.9, description is no longer parsed

from __future__ import print_function
import sys
import os
import screed


def main():
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} length-cutoff <seqfile> <longseqfile>' 
        print(mes.format(os.path.basename(sys.argv[0])), file=sys.stderr)
        sys.exit(1)

    min_length = int(float(sys.argv[1]))
    seqfile = sys.argv[2]
    longfile = sys.argv[3]
    if longfile == '-':
        longfile = '/dev/stdout'
    with open(longfile, 'w') as fw:
        for record in screed.open(seqfile):
            seq = record['sequence']
            if len(seq) >= min_length:
                print(
                    '>{}\n{}'.format(
                        record['name'], 
                        seq[:min_length]), 
                    file=fw)

if __name__ == '__main__':
    main()
