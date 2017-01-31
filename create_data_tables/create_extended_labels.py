#!/usr/bin/env python

"""
This script creates label tables from the csv label files.

Run it: 'python create_extended_labels.py TF1 TF2 ... TFN '

If no TF is given default is to run it on all the final submission TF-s.

Basically it only for easier matching between x-y values,
join is very slow in pandas.
"""

import sys,subprocess
import pandas as pd

def load_labels(t_factor,label_dir):
    """Load the labels for a transcription factor for all cell lines."""
    labels_fname = label_dir + t_factor + ".train.labels.tsv.gz"
    #load the original
    temp_df=pd.read_table(labels_fname,header=0,index_col=(0,1,2))
    #change U,A,B to U,A->0 B->1!!
    temp_df[(temp_df=='U') | (temp_df=='A')]=0
    temp_df[temp_df=='B']=1
    return temp_df


#TF list as argument or default
tf_list=sys.argv[1:]
if tf_list==[]:
    print 'No TF list given, running on all final submission TFs'
    tf_list=['ATF2','CTCF','E2F1','EGR1','FOXA1','FOXA2',
             'GABPA','HNF4A','JUND','MAX','NANOG','REST','TAF1']
    

#load index
print 'Reading the test regions ...',
sys.stdout.flush()
idx_df=pd.read_csv(
    '/encodeChallenge_Data/annotations/test_regions.blacklistfiltered.bed.gz',
    header=None,
    index_col=(0,1,2),
    sep='\t')
idx_df.index.rename(['chr','start','stop'],inplace=True)
print 'Done'

#create output dir
subprocess.call(['mkdir','extended_labels'])

for tf in tf_list:
    print 'Creating label table for ',tf,'...',
    sys.stdout.flush()
    #load labels
    labels=load_labels(tf,'/encodeChallenge_Data/labels/')   
    # join with the index
    labels=idx_df.merge(labels,how='left',left_index=True,right_index=True)
    print ' saving it ... ',
    sys.stdout.flush()
    labels.to_hdf('extended_labels/'+tf+'_labels.hdf','labels')
    print 'Done'