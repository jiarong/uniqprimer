from __future__ import print_function, unicode_literals

import sys
import os

def RC(seq):
    if isinstance(seq, bytes):
        seq = seq.decode('utf-8')
    table = {ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'}
    rc = seq.translate(table)
    return rc[::-1]

def main():
    if len(sys.argv) != 2:
        mes = '*** python {} "DNAString"'
        print(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    s = sys.argv[1]
    s = s.upper()
    if not set(s).issubset(set('ATCG')):
        mes = (
            '*** Illegal letters found, '
            'DNA string can only have ATCG..'
            )
        print(mes)
        sys.exit(1)
    print(RC(s))

if __name__ == '__main__':
    main()
