"""

Functions to create a submission for the DREAM encode challenge

"""

import pandas as pd
import numpy as np
import xgboost as xgb
import subprocess


""" 
Cell line definitions form here:
    https://www.synapse.org/#!Synapse:syn6131484/wiki/402034
    
"""

def get_cell_lines(tf,phase='final'):
    """Return the cell lines for a transcription factor."""
    #train cell lines
    tr={}
    tr['CTCF']=['A549','H1-hESC','HeLa-S3','HepG2','IMR-90','K562','MCF-7']
    tr['E2F1']=['GM12878','HeLa-S3']
    tr['EGR1']=['GM12878','H1-hESC','HCT116','MCF-7']
    tr['FOXA1']=['HepG2']
    tr['FOXA2']=['HepG2']
    tr['GABPA']=['GM12878','H1-hESC','HeLa-S3','HepG2','MCF-7']
    tr['HNF4A']=['HepG2']
    tr['JUND']=['HCT116', 'HeLa-S3', 'HepG2', 'K562', 'MCF-7']
    tr['MAX']=['A549','GM12878','H1-hESC','HCT116','HeLa-S3','HepG2','K562']
    tr['NANOG']=['H1-hESC']
    tr['REST']=['H1-hESC','HeLa-S3','HepG2','MCF-7','Panc1']
    tr['TAF1']=['GM12878', 'H1-hESC', 'HeLa-S3', 'K562']
    
    #test cell lines
    te={'final':{},'lb':{}}
    te['final']['CTCF']=['PC-3','induced_pluripotent_stem_cell']
    te['lb']['CTCF']=['GM12878']
    te['final']['E2F1']=['K562']
    te['lb']['E2F1']=None
    te['final']['EGR1']=['liver']
    te['lb']['EGR1']=['K562']
    te['final']['FOXA1']=['liver']
    te['lb']['FOXA1']=['MCF-7']
    te['final']['FOXA2']=['liver']
    te['lb']['FOXA2']=None
    te['final']['GABPA']=['liver']
    te['lb']['GABPA']=['K562']
    te['final']['HNF4A']=['liver']
    te['lb']['HNF4A']=None
    te['final']['JUND']=['liver']
    te['lb']['JUND']=['H1-hESC']
    te['final']['MAX']=['liver']
    te['lb']['MAX']=['MCF-7']
    te['final']['NANOG']=['induced_pluripotent_stem_cell']
    te['lb']['NANOG']=None
    te['final']['REST']=['liver']
    te['lb']['REST']=['K562']
    te['final']['TAF1']=['liver']
    te['lb']['TAF1']=['HepG2']
    
    return tr[tf],te[phase][tf]


""" Load data tables which were create before."""

def load_data(tf,phase):
    """
    Load all data for a transcription factor.
    
    It loads dnase fold change, deepbind score,
    baseline motifs, chipmunk/hocomoco/sarus motifs,
    and the labels.
    The loaded tables were created with the create_table...
    scripts.
    
    There are sequence realated tables, they are universal.
    The labels and the dnase tables are specific for a 
    cell line. The labels are only avaiable for the train
    cell lines. The dnase is available for train and test
    and therefore I am returning a list on dnase tables
    for the train cell lines and one table for test cell line.
    """
    #get training and test cell lines
    train_cl_list,test_cl_list=get_cell_lines(tf,phase)
    
    # load fold coverage tables for all train cell lines
    fc_train=[]
    for cl in train_cl_list:
        fn='/tables/dnase_fc/'+cl+'_dnase_fold_cov.hdf'
        fc_train.append(pd.read_hdf(fn))
    # load fold coverage tables for all test cell lines
    fc_test=[]
    for cl in test_cl_list:
        fn='/tables/dnase_fc/'+cl+'_dnase_fold_cov.hdf'
        fc_test.append(pd.read_hdf(fn))
    
    #load deepbind table
    deepbind=pd.read_hdf('/tables/deepbind/'+tf+'.hdf5')
    #load baseline motif table
    motif=pd.read_hdf('/tables/motif/'+tf+'_motif.hdf')
    
    #load sarus pwm motif table
    if tf!='TAF1':
        sarus=pd.read_hdf('/tables/sarus_pwm/'+tf+'_pwm.hdf')
    else:
        sarus = None

    #load labels for the transcription factor
    labels=pd.read_hdf('/tables/labels/'+tf+'_labels.hdf')
    
    #get a chromosome array for indexing
    #i had some serious performance issues with
    #pandas indexing on the 60 million row dataframes
    # so i am using numpy for the indexing later
    chr_arr=deepbind.reset_index()['chr'].values

    return fc_train,fc_test,deepbind,motif,sarus,labels,chr_arr


"""Stack the different inputs, labels into x,y data matrices."""

def create_tables(tf,train_chrs,valid_chrs,phase='final'):
    """
    Create data tables.
    
    Actually return dicionaries.They keys are the
    cell line names and the values are the numpy matrices.
    """
    #load cell lines
    train_cl_list,test_cl_list=get_cell_lines(tf,phase)
    
    #load all data
    fc_tr,fc_te,deepb,mot,sar,lab,ch_arr = load_data(tf,phase)
    
    #create train
    X_train,Y_train={},{}
    for cl,cl_fc in zip(train_cl_list,fc_tr):
        X_train[cl]=create_x(train_chrs,ch_arr,cl_fc,deepb,mot,sar)
        Y_train[cl]=create_y(train_chrs,ch_arr,cl,lab)
    #create valid
    X_valid,Y_valid={},{}
    for cl,cl_fc in zip(train_cl_list,fc_tr):
        X_valid[cl]=create_x(valid_chrs,ch_arr,cl_fc,deepb,mot,sar)
        Y_valid[cl]=create_y(valid_chrs,ch_arr,cl,lab)
    #create test
    X_test={}
    for cl,cl_fc in zip(test_cl_list,fc_te):
        if phase=='final':
            X_test[cl]=create_x(None,ch_arr,cl_fc,deepb,mot,sar)
        else:
            lb_chrs=['chr1','chr8','chr21']
            X_test[cl]=create_x(lb_chrs,ch_arr,cl_fc,deepb,mot,sar)
            
    return X_train,X_valid,X_test,Y_train,Y_valid

def create_x(chrs,chr_arr,fc,deepbind,motif,sarus):
    """Creat input data matrix."""
    # no restriction from chromosomes
    if chrs is None:
        cols=[fc.values,deepbind.values,motif.values]
        if sarus is not None:
            cols.append(sarus.values)
        x=np.column_stack(cols)
        return x
    else:
        ch_idx = get_chr_index(chr_arr,chrs)
        cols=[fc.values[ch_idx],
               deepbind.values[ch_idx],
               motif.values[ch_idx]]
        if sarus is not None:
            cols.append(sarus.values[ch_idx])
        x=np.column_stack(cols)
        return x

def create_y(chrs,chr_arr,cl,labels):
    """Creat labels data matrix."""
    ch_idx = get_chr_index(chr_arr,chrs)
    y=labels[cl].values[ch_idx].astype(int)
    return y

def get_chr_index(ch_arr,ch_list):
    """
    Create numpy bool chr mask.
    Pandas indexing has weird behaviour.
    """
    idx=np.full(len(ch_arr),False,dtype=bool)
    for ch in ch_list:
        idx|=ch_arr==ch
    return idx



""" 
Train xgboost ensemble.

The input of the model is not exactly the x created before
but a shifted version of x which contains the values for the 
neighborhood too. The shifted matrices are large and therefore
they are created on the fly.
"""

def train_xgb_ensemble(params,X_train,X_valid,X_test,
                       Y_train,Y_valid,window=3):
    """ Train xgb ensemble."""
    #train cell lines
    cl_tr_list=X_train.keys()
    # validation cell lines are rolled
    cl_val_list=cl_tr_list[1:]+[cl_tr_list[0]]
    
    y_pred={}
    for cl_tr,cl_val in zip(cl_tr_list,cl_val_list):
        print '\n-----------------------'
        print 'learning from '+cl_tr
        model=fit_xgb(params,X_train[cl_tr],Y_train[cl_tr],
                      X_valid[cl_val],Y_valid[cl_val],window)
        y_pred[cl_tr]=predict_xgb(model,X_test,window)
    return y_pred
    
def fit_xgb(params,x_train,y_train,x_valid,y_valid,window):
    """Fit xgb and predict."""
    #shift
    nx_train,ny_train=shift(x_train,y_train,window)
    nx_valid,ny_valid=shift(x_valid,y_valid,window)
    #create xgb matrices
    d_train = xgb.DMatrix(nx_train, label=ny_train)
    d_valid = xgb.DMatrix(nx_valid,label=ny_valid)
    #define printed evaluations
    evallist  = [(d_train,'train'),(d_valid,'eval')]
    #train
    bst = xgb.train(params,d_train,evals=evallist,
                    num_boost_round=200,
                    early_stopping_rounds=20,
                    verbose_eval=20)
    return bst

def predict_xgb(model,x,window):
    """Predict with xgb model."""
    #predict for all cell lines
    y_pred={}
    for cl_te,x_cl in x.iteritems():
        #predict
        y_pred[cl_te]=predict_on_batches_xgb(model,x_cl,window)
    return y_pred

def predict_on_batches_xgb(model,x,window,n=20):
    """Pedict xgb model on batches."""
    y_pred=[]
    chunk=int(1e6)
    
    #pad it for shift overlap in chunks borders
    if window!=1:
        x_pad=np.concatenate([x[-window/2:],x,x[:window/2]])
    else:
        x_pad=x
    
    for i in range(0,len(x_pad)-window/2,chunk-window/2):
        start,end=i,min(i+chunk+window/2,len(x_pad)-window/2)
        x_tmp=x_pad[start:end]
        if len(x_tmp)<window:
            continue
        if window!=1:
            x_tmp,_=shift(x_tmp,np.zeros(len(x_tmp)),window)
        d = xgb.DMatrix(x_tmp)
        y_pred.append(model.predict(d))
    y_pred = np.concatenate(y_pred)
    return y_pred

def shift(x,y,window):
    """Shift for neighborhood."""
    assert window%2==1
    x_l=[]
    for i in xrange(window):
        x_l.append(x[i:][:len(x)-window+1])
    nx=np.column_stack(x_l)
    ny=y[window/2:][:len(x)-window+1]
    return nx,ny


"""Ensembling, a simple average."""

def average_predictions(y_pred,tf,phase='final'):
    """Average predictions for all test cell lines."""
    _,test_cl_list=get_cell_lines(tf,phase)
    final_preds={}
    for test_cl in test_cl_list:
        y_all=[y[test_cl] for y in y_pred.values()]
        final_preds[test_cl]=np.mean(y_all,axis=0)
    return final_preds


"""Write the formatted submission."""

def write_subs(y_final,tf,phase='final'):
    """Write final submission."""
    #select annotation file and prefix for phase
    annot_dir='/encodeChallenge_Data/annotations/'
    if phase=='final':
        annot_file=annot_dir+'test_regions.blacklistfiltered.bed.gz'
        prefix='F.'
    else:
        annot_file=annot_dir+'ladder_regions.blacklistfiltered.bed.gz'
        prefix='L.'    
        
    #loop over tets cell lines
    _,test_cl_list=get_cell_lines(tf,phase)
    for test_cl in test_cl_list:
        write_sub(y_final[test_cl],prefix,tf,test_cl,annot_file)
    
def write_sub(y,prefix,tf,cl,annot_fn):
    """Write a submission"""
    #format the submission filename
    sub_fn=prefix+tf+'.'+cl+'.tab.gz'
    #save predicted values in tmp file
    tmp_fn=tf+'_'+cl+'.txt'
    with open(tmp_fn,'w') as f:
        f.write('\n'.join(y.astype('str')))
    #paste annotations
    cmd = ' paste   <( zcat ' + annot_fn +' ) '
    cmd+= tmp_fn + ' | gzip -c -1  > ' + sub_fn
    subprocess.check_output(cmd,shell=True,executable='/bin/bash')
    #rm tmp file
    subprocess.check_output(['rm',tmp_fn])
    