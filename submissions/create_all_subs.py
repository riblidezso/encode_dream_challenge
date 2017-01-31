#!/usr/bin/env python2

"""Create all submissions!!"""

import sys
import subprocess

phase=sys.argv[1]
window=sys.argv[2]

tf_l=dict()
tf_l['final']=['CTCF', 'E2F1', 'EGR1', 'FOXA1', 'FOXA2',
      'GABPA', 'HNF4A', 'JUND', 'MAX', 'NANOG',
      'REST', 'TAF1']
tf_l['lb']=['CTCF', 'EGR1', 'FOXA1', 'FOXA2',
      'GABPA', 'HNF4A', 'JUND', 'MAX', 'NANOG',
      'REST', 'TAF1']

for tf in tf_l[phase]:
    subprocess.call(['./create_submission.py',tf,phase,window])