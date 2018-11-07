#!/usr/bin/env python
from __future__ import print_function, unicode_literals

import sys
import os
import numpy
import yaml

import screed
import primer3

def RC(seq):
    if isinstance(seq, bytes):
        seq = seq.decode('utf-8')
    table = {ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'}
    rc = seq.translate(table)
    return rc[::-1]

def process(configd, desc_d, name, d_pos):

    l_f = sorted(set(d_pos['f']))
    l_r = sorted(set(d_pos['r']))
    a_f = numpy.array(l_f)
    a_r = numpy.array(l_r)
    range_st = range(configd['LEN_LOWER'], configd['LEN_UPPER'])
    for i, n in enumerate(a_f):
        f_name = '{}__{}'.format(name, n) 
        f_d = desc_d['f'][f_name]
        f = f_d['seq']
        for j, m in enumerate(l_r):
            diff = m - n
            if diff < 0:
                continue
            elif diff in range_st:
                r_name = '{}__{}'.format(name, m) 
                r_d = desc_d['r'][r_name]
                r = r_d['seq']
                # desc e.g: 
                #contiglen__184765;di__'f';tm__56.9135107847;gc__0.6
                if configd['SS']:
                    # check heterodimer
                    hetero = primer3.calcHeterodimer(
                        f.encode('utf-8'),
                        r.encode('utf-8'),
                        )
                    if hetero.dg < configd['DI_DG_LIMIT']:
                        continue
                # forward, f_tm, f_gc, reverse, r_tm, r_gc
                mes = (
                    '{}\t{}\t'
                    '{}\t{}\t{}\t'
                    '{}\t{:.1f}\t{:.2f}\t'
                    '{}\t{:.1f}\t{:.2f}\t'
                    )
                print(
                    mes.format(
                        name, f_d['contiglen'], 
                        n, m, diff,
                        f, float(f_d['tm']), float(f_d['gc']), 
                        r, float(r_d['tm']), float(r_d['gc']),
                        )
                    )

            # break when out of range since already sorted
            elif diff > configd['LEN_UPPER']:
                break


def main():
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} params.config file.filtindi'
        print(
            mes.format(os.path.basename(sys.argv[0])), 
            file=sys.stderr
            )
        sys.exit(1)

    configf = sys.argv[1]
    infile = sys.argv[2]
    configd = yaml.load(open(configf))


    print(
        '#contig_name\tcontig_len\t'
        'f_start\tr_start\tamp_size\t'
        'f_seq\tf_tm\tf_gc\t' 
        'r_seq\tr_tm\tr_gc\t'
        )
    for rec in screed.open(infile):
        name, desc = rec.name.split(None, 1)  #1095362_contig_14__2939
        # desc e.g: contiglen__184765;tm__56.9135107847;gc__0.6

        #dict(i.split('__') for i in desc.split(';'))
        _l = name.split('__')
        if len(_l) < 3:
            print(
                '*** {} does have wrong format, skipping..'.format(
                    name
                    ), 
                file=sys.stderr,
                )
            continue
        contig_name = _l[0]
        pos = int(_l[1])
        di = _l[2] # sense: 'f' ; anti-sense: 'r'
        name = '{}__{}'.format(contig_name, pos)
        try:
            if contig_name == prev_contig_name:
                d_pos[di].append(pos)
                d[di][name] = dict(
                    i.split('__') for i in desc.split(';')
                    )
                d[di][name]['seq'] = rec.sequence

            else:
                process(configd, d, prev_contig_name, d_pos)
                # initialization 
                del d_pos
                del d
                d = {'f':{}, 'r':{}}
                d[di][name] = dict(
                    i.split('__') for i in desc.split(';')
                    )
                d[di][name]['seq'] = rec.sequence
                d_pos = {'f':[], 'r':[]}
                d_pos[di].append(pos)
                prev_contig_name = contig_name

        except UnboundLocalError as e:
            # for first record
            #UnboundLocalError due to prev_contig_name not defined
            d = {'f':{}, 'r':{}}
            d[di][name] = dict(
                i.split('__') for i in desc.split(';')
                )
            d[di][name]['seq'] = rec.sequence
            d_pos = {'f':[], 'r':[]}
            d_pos[di].append(pos)
            prev_contig_name = contig_name
    # process last batch
    try:
        process(configd, d, prev_contig_name, d_pos)
    except UnboundLocalError as e:
        mes = (
            '*** Empty file detected: {} '
            '(file.uniq2ref.primer.filterindi), skipping..'
            )
        print(
            mes.format(os.path.basename(infile)),
            file=sys.stderr,
            )
        sys.exit(0)

if __name__ == '__main__':
    main()
