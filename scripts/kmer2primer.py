#!/usr/bin/env python
from __future__ import print_function

import sys
import os

import yaml

import screed

def parse_kmerfile(f):
    d = {}
    for rec in screed.open(f):
        name = rec.name.split(None, 1)[0]  #1095362_contig_14__2939
        _l = name.split('__')
        if len(_l) < 2:
            print('*** {} does have wrong format, skipping..'.format(name), 
                    file=sys.stderr)
            continue
        contig_name = _l[0]
        pos = int(_l[1])
        _st = d.get(contig_name, set())
        d[contig_name] = _st
        _st.add(pos)  # this way due to shallow copy

    return d

def main():
    if len(sys.argv) != 4:
        mes = (
            '*** Usage: python {} params.config '
            'concontig.fasta file.uniq2ref'
            )
        print(
            mes.format(os.path.basename(sys.argv[0])), 
            file=sys.stderr,
            )
        sys.exit(1)

    configf = sys.argv[1]
    contigf = sys.argv[2]
    kmerfile = sys.argv[3]

    configd = yaml.load(open(configf))
    d = parse_kmerfile(kmerfile)
    if len(d) == 0:
        mes = '*** Empty file detected: {} (file.uniq2ref), skipping..'
        print(
            mes.format(os.path.basename(kmerfile)),
            file=sys.stderr,
            )
        sys.exit(0)

    for rec in screed.open(contigf):
        name = rec.name
        name = name.split(None, 1)[0]
        if not name in d:
            continue
        seq = rec.sequence

        for p in d[name]:
            assert len(seq) > p, (
                '*** seq length < primer start position'
                )
            primer = seq[p : p+configd['K']]
            _mes = '>{}__{}  contig__{}\n{}'
            print(_mes.format(name, p, len(seq), primer), file=sys.stdout)

if __name__ == '__main__':
    main()
