#!/usr/bin/env python

"""
This script creates DNASE fold covarage tables from the bigwig files.

Run it: 'python create_fold_cov_table.py CELL_LINE_1 CELL_LINE_2 ... CELL_LINE_N '

If no cell line is given default is to run it on all the cell lines.
"""

import os,subprocess,sys,time
from multiprocessing import Pool
from multiprocessing import sharedctypes
from numpy import ctypeslib
from functools import partial
import numpy as np
import pandas as pd
import pyBigWig

DNASE_FOLD_COV_DIR ="/encodeChallenge_Data/DNASE/fold_coverage_wiggles/"

def load_dnase_fold_cov(cell_line,
                        n_proc=12,
                        dnase_fold_cov_dir=DNASE_FOLD_COV_DIR ):
    """Load the raw fold coverage for a t factor and cell line."""
    dnase_fold_cov_fname = dnase_fold_cov_dir + 'DNASE.'+cell_line+'.fc.signal.bigwig'
    
    #partially apply the scorer function
    part_get_dnase_fold_cov_from_region=partial(get_dnase_fold_cov_from_region,
                                                fn=dnase_fold_cov_fname)
    #parallel execute it
    Pool(n_proc).map(part_get_dnase_fold_cov_from_region,xrange(len(idx)))

    #return a dataframe
    return pd.DataFrame({cell_line+'_dnase_fc' : fold_cov},index=idx)

def get_dnase_fold_cov_from_region(i,fn=None):
    """Get raw dnase fold coverage for a region."""
    contig,start,stop=idx[i]
    bigwig_f=pyBigWig.open(fn)
    
    fc=bigwig_f.stats(contig,start,stop)[0]
    if fc!=None:
        fold_cov[i]=fc
    else:
        fold_cov[i]=0
        
    bigwig_f.close()
    return



###########################

#argv is the cell line list
cl_list=sys.argv[1:]
if cl_list==[]:
    print 'No cell line list given, running on all cell lines'
    cl_list=[x.split('.')[1] for x in os.listdir(DNASE_FOLD_COV_DIR)]

#load index
print 'Reading the test regions ...',
sys.stdout.flush()
idx=pd.read_csv(
    '/encodeChallenge_Data/annotations/test_regions.blacklistfiltered.bed.gz',
    header=None,
    index_col=(0,1,2),
    sep='\t').index
idx.rename(['chr','start','stop'],inplace=True)
print 'Done'

#global shared result array... is there a better way?
fold_cov = sharedctypes.RawArray('d', len(idx))

#create output dir
subprocess.call(['mkdir','fold_cov_data'])

for cl in cl_list:  
    print 'Creating fold cov table for ',cl,'...',
    dnase_fc_df=load_dnase_fold_cov(cl)
    print ' saving it ... ',
    dnase_fc_df.to_hdf('fold_cov_data/'+cl+'_dnase_fold_cov.hdf',
                       'dnase_fold_cov')
    print 'Done'
