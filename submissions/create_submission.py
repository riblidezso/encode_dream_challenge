#!/usr/bin/env python2

"""
Command line tool to create a submission.

Both final and leaderboard submissions can be created.
"""

import sys
import time
import submission_utils as sut

#benchmark exec time
start=time.time()

#which submission
tf= sys.argv[1]
phase= sys.argv[2]
#size of neighborhood
#recommended 3?
window = int(sys.argv[3])

#chromosomes used for training and validation
#with early stopping
chr_tr= ['chr3','chr4','chr6','chr7','chr19','chr20']
#['chr3','chr4','chr5','chr6','chr7','chr10','chr11',
#        'chr12','chr13','chr14','chr15','chr16','chr17',
#        'chr18','chr19','chr20']
chr_va=['chr9','chr22']
#['chr2','chr9','chr22']

#xgboost model params
params = {'max_depth':5,
          'eta':0.1,
          'min_child_weight':10,
          'colsample_bytree':0.7,
          'subsample':0.8,
          'silent':1,
          'objective': "binary:logistic",
          'eval_metric': 'auc',
          'nthread':20}

#load data
X_tr,X_va,X_te,Y_tr,Y_va= sut.create_tables(tf,chr_tr,chr_va,phase)

#train xgb ensemble
y_pred=sut.train_xgb_ensemble(params,X_tr,X_va,X_te,Y_tr,Y_va,window)

#average predictions
y_final=sut.average_predictions(y_pred,tf,phase)

# write submission file
sut.write_subs(y_final,tf,phase)

print 'Done:',int((time.time()-start)/60),'minutes''chr18'