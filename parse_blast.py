#!/usr/bin/env python3.4

import sys
seqid = ''

with open(sys.argv[1]) as f:
    for line in f:
        if line.split('\t')[0] != seqid:
            print(line, end='')
            seqid = line.split('\t')[0]
