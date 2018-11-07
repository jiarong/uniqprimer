#!/usr/bin/env python
from __future__ import print_function

import sys
import os

def parse(fp):
    pair = False
    for line in fp:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('# EPRIMER32 RESULTS FOR'):
            name = line.split()[-1]
            continue
        elif line.startswith('#'):
            continue

        line = line.strip()
        if 'PRODUCT SIZE:' in line:
            pair = True
            size = int(line.split()[-1])
            continue

        elif 'FORWARD PRIMER' in line:
            #FORWARD PRIMER  542890   20  59.97  55.00  CGCGCTTGAATAGTCGTTGG
            _, _, s_f, Len_f, Tm_f, GC_f, Seq_f = line.split()
            assert pair == True, '*** Check format..'

            continue

        elif 'REVERSE PRIMER' in line:
            #FORWARD PRIMER  542890   20  59.97  55.00  CGCGCTTGAATAGTCGTTGG
            _, _, s_r, Len_r, Tm_r, GC_r, Seq_r = line.split()
            assert pair == True, '*** Check format..'
            assert int(s_r) - int(s_f) + int(Len_r) == size
            mes = ('>{}__{}to{}/1  Len_{};Tm_{};GC_{};size_{}\n{}\n'
                    '>{}__{}to{}/2  Len_{};Tm_{};GC_{};size_{}\n{}')
            print(mes.format(name, s_f, s_r, Len_f, Tm_f, GC_f, size, Seq_f,
                                name, s_f, s_r, Len_r, Tm_r, GC_r, size,
                                Seq_r))
            pair = False

if __name__ == '__main__':
    if len(sys.argv) != 2:
        mes = '*** Usage: python {} file.ebeprimer'
        print(mes.format(os.path.basename(sys.argv[0])), file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    try:
        if infile == '-':
            fp = sys.stdin
        else:
            fp = open(infile)

        parse(fp)
    finally:
        fp.close()
