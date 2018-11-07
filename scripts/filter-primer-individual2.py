#!/usr/bin/env python
from __future__ import print_function, unicode_literals

import sys
import os

import yaml

import screed
import primer3
import numpy

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
    cnt = 0
    for rec in screed.open(primerfile):
        cnt += 1
        _name = rec.name
        name, _contig = _name.split(None, 1)
        contig_len = _contig.split('__', 1)[1]
        seq = rec.sequence
        # primer3 functions only accept byte-strings
        seq = seq.encode('utf-8')
        #seq = bytes(seq, 'utf-8')
        seq_rc = RC(seq)

        a_ambi = numpy.array(has_ambiguous(seq), has_ambiguous(seq_rc))
        if sum(a_ambi) == 2:
            continue 

        # check tm
        tm = primer3.calcTm(seq)
        tm_rc = primer3.calcTm(seq_rc)
        a_tm = numpy.array(
            (tm < d['TM_LOWER'] or tm > d['TM_UPPER']),
            (tm_rc < d['TM_LOWER'] or tm_rc > d['TM_UPPER']),
            )
        if sum(a_tm) == 2:
            continue

        # check gc
        gc = check_gc(seq)
        gc_rc = check_gc(seq_rc)
        a_gc = numpy.array(
            (gc < d['GC_LOWER'] or gc > d['GC_UPPER']),
            (gc_rc < d['GC_LOWER'] or gc_rc > d['GC_UPPER']),
            )
        if sum(a_gc) == 2:
            continue

        if d['GC_CLAMP']:   
            c = end_gc_count(seq)
            c_rc = end_gc_count(seq_rc)
            a_endgc = numpy.array(
                c > 3 or c < 1,
                c_rc > 3 or c_rc < 1,
                )
            if sum(a_endgc) == 2:
                continue

        if d['SS']:
            hp = primer3.calcHairpin(seq)
            ho = primer3.calcHomodimer(seq)
            hp_rc = primer3.calcHairpin(seq_rc)
            ho_rc = primer3.calcHomodimer(seq_rc)
            orig_pass = (
                (hp.dg < d['HP_DG_LIMIT'] or hp.dg > 0)
                & (ho.dg < d['DI_DG_LIMIT'] or ho.dg > 0)
                )
            rc_pass = (
                (hp_rc.dg < d['HP_DG_LIMIT'] or hp_rc.dg > 0)
                & (ho_rc.dg < d['DI_DG_LIMIT'] or ho_rc.dg > 0)
                )
            if ho.dg < d['DI_DG_LIMIT'] or ho.dg > 0:
                continue

        pass_cnt += 1
        mes = '>{}  contiglen__{};tm__{};gc__{}\n{}'
        print(mes.format(name, contig_len, tm, gc, seq), file=sys.stdout)

    if cnt == 0:
        mes = '*** Empty file detected: {} (file.uniq2ref.primer), skipping..'
        print(
            mes.format(os.path.basename(primerfile)),
            file=sys.stderr,
            )
        sys.exit(0)

if __name__ == '__main__':
    main()
