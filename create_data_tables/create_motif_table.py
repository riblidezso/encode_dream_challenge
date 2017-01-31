#!/usr/bin/env python

"""
This script creates TF motif tables.

Run it: 'python create_motif_table.py TF1 TF2 TF3 ... TFN '

If no TF is given default is to run it on all the final submission round TF-s.
"""


import os,subprocess,sys,time
from multiprocessing import Pool,sharedctypes
from numpy import ctypeslib
from functools import partial
import numpy as np
import pandas as pd
from scipy.stats.mstats import mquantiles
from Bio import SeqIO
from pyDNAbinding.binding_model import DNASequence, PWMBindingModel, DNABindingModels, load_binding_models



def aggregate_region_scores(scores,
                            quantile_probs = [0.99, 0.95, 0.90, 0.75, 0.50]):
    """Return aggregate scores of all scores from a region."""
    rv = [scores.mean()/len(scores), scores.max()]
    rv.extend(mquantiles(scores, prob=quantile_probs))
    return rv

def load_motif_scores(t_factor,
                      n_proc=12):
    """Load aggregate motif scores for a t_factor for intervals in the index."""
    #load binding models
    binding_models = load_binding_models("models.yaml")
    model = binding_models.get_from_tfname(t_factor)
    
    #partially apply the scorer function
    part_get_motif_scores_from_region=partial(
        get_motif_scores_from_region,model=model)
    
    #parrallel execute it
    Pool(n_proc).map(part_get_motif_scores_from_region,xrange(len(idx)))        
    
    #shared array to numpy array
    ms=ctypeslib.as_array(motif_scores).reshape(shape)
    
    colnames=[t_factor+'_'+x for x in aggregate_region_scores_labels]
    return pd.DataFrame(ms,index=idx,columns=colnames)
    
def get_motif_scores_from_region(i,model):
    """Get motif scores from a region in index."""
    contig, start, stop = idx[i]
    #load seq from genome
    seq=str(ref_genome[contig][start:stop].seq).upper()
    
    #make sharred array numpy array
    ms = ctypeslib.as_array(motif_scores).reshape(shape)
    
    #get aggregate motif scores
    ms[i] = aggregate_region_scores(
        DNASequence(seq).score_binding_sites(model[0], 'MAX'))
    
    return


###########################

#TF list as argument or default
tf_list=sys.argv[1:]
if tf_list==[]:
    print 'No TF list given, running on all final submission TFs'
    tf_list=['ATF2','CTCF','E2F1','EGR1','FOXA1','FOXA2',
             'GABPA','HNF4A','JUND','MAX','NANOG','REST','TAF1']

#load ref genome
print 'Reading reference genome ...',
sys.stdout.flush()
REF_GENOME="/encodeChallenge_Data/hg19.genome.fa"
with open(REF_GENOME, "rU") as h:
    ref_genome=SeqIO.to_dict(SeqIO.parse(h, "fasta"))
print 'Done'

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


#global var for baselines script
#the aggregate colnames
aggregate_region_scores_labels = ["motif_mean", "motif_max", "motif_q99",
                                  "motif_q95", "motif_q90", "motif_q75", "motif_q50"]

#global var for my multiprocessing
#global shared result array... is there a better way?
shape=(len(idx),len(aggregate_region_scores_labels))
motif_scores = sharedctypes.RawArray('d', shape[0]*shape[1])

#create output dir
subprocess.call(['mkdir','motif_data'])

#loop over the TFs
for tf in tf_list:
    #create motif score table
    print 'Creating motif table for ',tf,'...',
    sys.stdout.flush()
    motif_df=load_motif_scores(tf)
    #rename index for
    #save it as hdf
    print ' saving it ... ',
    sys.stdout.flush()
    motif_df.to_hdf('motif_data/'+tf+'_motif.hdf','motif')
    print 'Done'