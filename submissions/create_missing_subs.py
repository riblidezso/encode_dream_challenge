#!/usr/bin/env python2

"""Create all submissions!!"""

import sys
import subprocess

phase=sys.argv[1]
window=sys.argv[2]

tf_l=dict()
tf_l['final']=['E2F1','FOXA2']

for tf in tf_l[phase]:
    subprocess.call(['./create_submission.py',tf,phase,window])