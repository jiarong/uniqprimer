#!/usr/bin/env python
from __future__ import print_function, unicode_literals

import os
import sys

import screed
import primer3

# python <thisfile> <check-kmer-distance-py-output> <contigs.fa>

# Refs for parameters:
#   http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
#   http://www.genomecompiler.com/tips-for-efficient-primer-design/
#   https://www.idtdna.com/pages/support/faqs/how-do-you-calculate-the-annealing-temperature-for-pcr-

TM_LOWER=57
TM_UPPER=62
TM_DIFF_MAX=5

GC_LOWER=0.4
GC_UPPER=0.6

SS=True  # secondary structures
HP_DG_LIMIT=-4000  # hairpin: 2K at 3' and 3K at mid (cal/mol)
DI_DG_LIMIT=-7000  # dimers 5K at 3' and 6K at mid

def parse_kmer_distance(f):
    try:
        if f == '-':
            fp = sys.stdin
        else:
            fp = open(f)
        d = {}
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line:
                continue
            name, f, r = line.split()
            try:
                d[name].add((int(f), int(r)))
            except KeyError as e:
                d[name] = set()
                d[name].add((int(f), int(r)))

        return d

    finally:
        fp.close()

def check_GC(seq):
    total = len(seq)
    count = 0
    for c in seq:
        if c in 'GC':
            count += 1
    return count*1.0/total

def has_ambiguous(seq):
    c_st = set(seq)
    if c_st.issubset(set('ATCG')):
        return False
    else:
        return True

def RC(seq):
    rc = seq.translate({ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'})
    return rc[::-1]
    
def main():
    if len(sys.argv) != 4:
        mes = '*** python {} size <check-kmer-distance-py-output> <contigs.fa>'
        print(mes.format(os.path.basename(sys.argv[0])), file=sys.stderr)
        sys.exit(1)

    size = int(sys.argv[1])
    infile = sys.argv[2]
    contigf = sys.argv[3]

    d = parse_kmer_distance(infile)

    print(('#contig_name\tcontig_len\tf_start\tr_start\tf_seq\tf_tm\tf_gc\t'
             'r_seq\tr_tm\tr_gc\tta\tamp_size'))
    pair_pass = 0
    for rec in screed.open(contigf):
        name = rec.name
        if not name in d:
            continue
        seq = rec.sequence

        for f_p, r_p in d[name]:
            assert len(seq) > r_p, '*** seq length < forward primer position'
            f = seq[f_p:(f_p+size)]
            r = RC(seq[r_p:(r_p+size)])

            # primer3 functions only accept byte-strings
            f = f.encode('utf-8')
            #f = bytes(f, 'utf-8')
            r = r.encode('utf-8')
            #r = bytes(r, 'utf-8')
            if has_ambiguous(f) or has_ambiguous(r):
                continue 

            # check tm
            f_tm = primer3.calcTm(f)
            if f_tm < TM_LOWER or f_tm > TM_UPPER:
                continue
            r_tm = primer3.calcTm(r)
            if r_tm < TM_LOWER or r_tm > TM_UPPER:
                continue
            if abs(f_tm - r_tm) > TM_DIFF_MAX:
                continue

            # check gc
            f_gc = check_GC(f)
            if f_gc < GC_LOWER or f_gc > GC_UPPER:
                continue
            r_gc = check_GC(r)
            if r_gc < GC_LOWER or r_gc > GC_UPPER:
                continue

            amp = seq[f_p:(r_p+size)].encode('utf-8')
            amp_tm = primer3.calcTm(amp)
            ta = 0.3*min(f_tm,r_tm) + 0.7*amp_tm - 14.9 # premierbiosoft
            #ta = 0.3*min(f_tm,r_tm) + 0.7*amp_tm - 25 # IDT recommendation
            ### thermodynamics check
            ### skipping here as loose filter
            # check hairpin and homodimer

            if SS:
                f_hp = primer3.calcHairpin(f)
                f_ho = primer3.calcHomodimer(f)


                if f_hp.dg < HP_DG_LIMIT or f_hp.dg > 0:
                    continue
                if f_hp.tm > ta:
                    continue 
                if f_ho.dg < DI_DG_LIMIT or f_ho.dg > 0:
                    #print('+++++>', f_ho.dg)
                    continue


                r_hp = primer3.calcHairpin(r)
                r_ho = primer3.calcHomodimer(r)
                if r_hp.dg < HP_DG_LIMIT or r_ho.dg > 0:
                    continue
                if r_hp.tm > ta:
                    continue 
                if r_ho.dg < DI_DG_LIMIT or r_ho.dg > 0:
                    #print('=====>', r_ho.dg)
                    continue

                # check heterodimer
                hetero = primer3.calcHeterodimer(f,r)
                if hetero.dg < DI_DG_LIMIT:
                    continue

            pair_pass += 1
            # forward, f_tm, f_gc, reverse, r_tm, r_gc, ta, amp_size
            mes = ('{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.2f}\t{}\t{:.1f}\t'
                     '{:.2f}\t{}\t{}')
            print(mes.format(name, len(seq), f_p, r_p, f, f_tm, f_gc, r, 
                              r_tm, r_gc, ta, len(amp)))

    print('*** Pairs passed: {}'.format(pair_pass), file=sys.stderr)

if __name__ == '__main__':
    main()
