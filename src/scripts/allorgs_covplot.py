#!/usr/bin/env python
"""This script is used to run covplot on all directories ending in _S* and for all organisms
present in orgs_list.csv (except those with parenthesis in the name, to avoid )
"""

import glob
import json
import os
import shlex
import subprocess
import sys

import pandas as pd

info = {}
oh = open("dump.json", "w")

for d in glob.glob("*_S?"):
    print(d)
    os.chdir(d)
    info[d] = {}
    orgs = pd.read_csv("orgs_list.tsv", delimiter="\t", header=0)
    scinames = set(orgs.ssciname.tolist())
    print(scinames)
    for sci in scinames:
        print(d, sci, file=sys.stderr)
        info[d][sci] = {}
        if "phage" in sci or "Phage" in sci:
            print("skipping %s" % sci)
            continue
        print("running covplot on %s" % sci, file=sys.stderr)
        cml = "virmet covplot --outdir ./ --organism '%s'" % sci.split("(")[
            0
        ].rstrip(" ")
        res = subprocess.check_output(shlex.split(cml))
        print(res)
        for keyval in res.decode("ascii").strip().split():
            print(keyval)
            key, val = str(keyval).split(":")
            info[d][sci][key] = val
    os.chdir("..")
json.dump(info, oh)
