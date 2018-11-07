from __future__ import print_function, unicode_literals

import sys
import os

import screed

def RC(seq):
    if isinstance(seq, bytes):
        seq = seq.decode('utf-8')
    table = {ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'}
    rc = seq.translate(table)
    return rc[::-1]

def parse_primer_summary(f):
    if f == '-':
        f = '/dev/stdin'
    with open(f) as fp:
        d = {}
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            contig_name = line.split(None,1)[0]
            try:
                d[contig_name].add(line)
            except KeyError as e:
                d[contig_name] = set([line,])

        return d

def main():
    if len(sys.argv) != 3:
        mes = (
            '*** python {} contig.fasta '
            'file.uniq2ref.primer.filtindi.filtpair'
            )
        print(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    contigf = sys.argv[1]
    sumfile = sys.argv[2]

    d = parse_primer_summary(sumfile)
    for rec in screed.open(contigf):
        name = rec.name.split(None, 1)[0]
        if not name in d:
            continue
        seq = rec.sequence
        for line in d[name]:
            _l =  line.split()
            assert len(_l) == 11, _l
            (_, contig_len, f_start, r_start, 
                amp_size, f_seq, f_tm, f_gc, r_seq, 
                r_tm, r_gc,) = _l

            assert seq.find(f_seq) == int(f_start)
            assert seq.find(RC(r_seq)) == int(r_start)
            assert len(seq) == int(contig_len)

    print('*** PASS =================', file=sys.stderr)

if __name__ == '__main__':
    main()
