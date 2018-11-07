#!/usr/bin/env python
from __future__ import print_function, unicode_literals

import sys
import os

import yaml

import screed
import primer3

def check_gc(seq):
    seq = seq.upper()
    total = len(seq)
    count = 0
    for c in seq:
        if c in 'GC':
            count += 1
    return count*1.0/total

def end_gc_count(seq):
    seq = seq.upper()
    count = 0
    for c in seq[-5:]:
        if c in 'GC':
            count += 1
    return count

def has_ambiguous(seq):
    c_st = set(seq)
    if c_st.issubset(set('ATCG')):
        return False
    else:
        return True

def RC(seq):
    if isinstance(seq, bytes):
        seq = seq.decode('utf-8')
    rc = seq.translate(
        {ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'}
        )
    return rc[::-1]

def main():
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} params.config file.uniq2ref.primer'
        print(
            mes.format(os.path.basename(sys.argv[0])), 
            file=sys.stderr,
            )
        sys.exit(1)

    configf = sys.argv[1]
    primerfile = sys.argv[2]

    d = yaml.load(open(configf))

    pass_cnt = 0
    total_cnt = 0
    for rec in screed.open(primerfile):
        total_cnt += 1
        _name = rec.name
        name, _contig = _name.split(None, 1)
        contig_len = _contig.split('__', 1)[1]
        seq = rec.sequence
        seq_rc = RC(seq)
        # primer3 functions only accept byte-strings
        seq = seq.encode('utf-8')
        seq_rc = seq_rc.encode('utf-8')
        #seq = bytes(seq, 'utf-8')
        trig = False
        for di, seq in zip(('f', 'r'), (seq, seq_rc)):
            if has_ambiguous(seq):
                continue 

            # check tm
            tm = primer3.calcTm(seq)
            if tm < d['TM_LOWER'] or tm > d['TM_UPPER']:
                continue

            # check gc
            gc = check_gc(seq)
            if gc < d['GC_LOWER'] or gc > d['GC_UPPER']:
                continue

            if d['GC_CLAMP']:   
                cnt = end_gc_count(seq)
                if cnt > 3 or cnt < 1:
                    continue

            if d['SS']:
                hp = primer3.calcHairpin(seq)
                ho = primer3.calcHomodimer(seq)
                if hp.dg < d['HP_DG_LIMIT'] or hp.dg > 0:
                    continue
                if ho.dg < d['DI_DG_LIMIT'] or ho.dg > 0:
                    continue

            trig = True
            mes = '>{}__{}  contiglen__{};di__{};tm__{};gc__{}\n{}'
            print(
                mes.format(name, di, contig_len, di, tm, gc, seq), 
                file=sys.stdout,
                )

        if trig:
            pass_cnt += 1

    mes = '*** # of primers (at least one direction) passed filter: {}' 
    print(mes.format(pass_cnt), file=sys.stderr)
    if total_cnt == 0:
        mes = (
            '*** Empty file detected: {} (file.uniq2ref.primer), '
            'skipping..'
            )
        print(
            mes.format(os.path.basename(primerfile)),
            file=sys.stderr,
            )
        sys.exit(0)

if __name__ == '__main__':
    main()
