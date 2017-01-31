#!/usr/bin/env python2

"""
Score a Transcription factor using deepbind.

Mulitprocessing added to speed up things.
Results are stored in [TF_NAME].hdf5 ad pd datafram in hdf5 format

execute it:
./create_deepbind.py [TF_NAME]
"""

N_PROC=20
REF_FASTA="/encodeChallenge_Data/hg19.genome.fa"
TEST_REGIONS='/encodeChallenge_Data/annotations/test_regions.blacklistfiltered.bed.gz'

import time,sys,os,subprocess
import random
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
import multiprocessing as mp

def predict_parall(tf,ref_genome,regions,n_split):
    """Predict region with deepbind in parallel."""
    #write id file
    id_fn=write_id_file(tf)
    
    #write the sequence files
    region_l=chunkIt(regions,n_split)
    seq_fns=[write_seq_f(ref_genome,r) for r in region_l]
    
    #define params
    res_fns=['dbres_'+str(i).zfill(3) for i in xrange(n_split)]
    params=zip(seq_fns,n_split*[id_fn],res_fns)
    
    #run it
    my_parall_runnner(predict_worker,params)
    
    #save results as hfd5
    save_hdf5(res_fns,regions,tf+'.hdf5')
    
    #remove tmp files
    for f1,f2 in zip(seq_fns,res_fns):
        subprocess.check_output(['rm',f1,f2])
    subprocess.check_output(['rm',id_fn])
    
     
def write_id_file(tf):
    """Wrtie id file for TF."""
    #load id db
    db=pd.read_csv('/common/tools/deepbind/db/db.tsv',
                   sep='\t',comment='#')
    #find ids
    ids=db[(db['Protein']==tf) &
       (db['Labels']!='deprecated')]['ID']
    #write id file
    id_fn=tf+'.ids'
    with open(id_fn,'w') as f:
        for idi in ids:
            f.write(idi+'\n')
    return id_fn

def chunkIt(seq, num):
    """Chunk a list into num parts"""
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def write_seq_f(ref_genome,regions):
    """Write a sequence file and return its name."""
    tmp_fn='tmp_'+str(random.randint(0,1e16))+'.seq'
    with open(tmp_fn,'w') as f:
        for ch,start,stop in regions:
            s=str(ref_genome[ch][start:stop].seq)
            f.write(s+'\n')
    return tmp_fn

def my_parall_runnner(f,args):
    """Run a function in parallel.""" 
    print 'Starting processes'
    procs=[mp.Process(target=f,args=arg) for arg in args]
    for proc in procs:
        proc.start()
        
    print '\nAll blocks started, waiting for the last ones to finish'
    sys.stdout.flush()
    finished = False
    while(not finished):
        finished = True
        for proc in procs:
            if (proc.is_alive() == True):
                finished = False
        time.sleep(0.1)

def predict_worker(seq_fn,id_fn,res_fn):
    """Predict deppbind on sequences.""" 
    cmd='/common/tools/deepbind/deepbind '
    cmd+=id_fn+' '+seq_fn+' > '+res_fn
    res=subprocess.call(cmd,shell=True)
    
def save_hdf5(res_fns,regions,out_fn):
    """Load results and save them as hdf5."""
    df_l=[pd.read_csv(fn,sep='\t') for fn in res_fns ]
    res=pd.concat(df_l).reset_index(drop=True)
    res=pd.concat([res,pd.DataFrame(regions)],axis=1)
    res.rename(columns={0:'chr',1:'start',2:'stop'},inplace=True)
    res.set_index(['chr','start','stop'],inplace=True)
    res.to_hdf(out_fn,'deepbind')


#####################################
if __name__ == "__main__":
    #get tf
    TF=sys.argv[1]
    
    print 'Reading reference genome ...',
    sys.stdout.flush()
    with open(REF_FASTA, "rU") as h:
        ref_genome=SeqIO.to_dict(SeqIO.parse(h, "fasta"))
    print
        
    print 'Reading the test regions ...',
    sys.stdout.flush()
    idx_df=pd.read_csv(
        TEST_REGIONS,
        header=None,
        sep='\t')
    regions=idx_df.values
    print

    print 'Scoring ...',
    start=time.time()
    predict_parall(TF,ref_genome,regions,N_PROC)
    print time.time()-start,'s'