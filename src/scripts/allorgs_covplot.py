#!/usr/bin/env python

import glob
import os
import pandas as pd
import subprocess
import shlex
import json

info = {}
oh = open('dump.json', 'w')

for d in glob.glob('*_S*'):
    print(d)
    os.chdir(d)
    info[d] = {}
    orgs = pd.read_csv('orgs_list.tsv', delimiter='\t', header=0)
    scinames = set(orgs.ssciname.tolist())
    print(scinames)
    for sci in scinames:
        info[d][sci] = {}
        if "(" in sci:
            print('skipping %s' % sci)
            continue
        print('running covplot on %s' % sci)
        cml = "virmet covplot --outdir ./ --organism '%s'" % sci
        res = subprocess.check_output(shlex.split(cml))
        for keyval in res.decode('ascii').split():
            print(keyval)
            key, val = str(keyval).split(':')
            info[d][sci][key] = val
        print(res)
    os.chdir('..')
    break
json.dump(info, oh)
