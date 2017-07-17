# Eotvos Lorand University team, ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge

Authors: Dezso Ribli, Janos Szalai-Gindl, Akos Rudas, Sandor Laki, Attila Kiss, Istvan Csabai

Affiliations: Eotvos Lorand University

The submission can be made public.

## Background/Intro


This is the first time we work with a transcription factor binding problem. We are physicist and computer scientists working with machine learning, data minin and various other subjects. We have tried to create a meaningful table representation of the data, using multiple transciption binding prediction methods, and apply machine learning on the table.


## Methods

#### Features

We have used the following variables:
- Motif scores as calculated in the baseline method. Used parallelism to speed up things.
- Dnase raw fold coverage values for every cell line in the position. Also parallelized.
- Deepbind scores.
- Sarus score from the hocomoco database.
- We use 3 region length sliding window to incorporate data from the neighbourhood.
- The predicted labels were transformed as the following (U->0,A->0,B->1)


Please note that our leaderboard submission is different from the final submissions, as it does not use the chipseq labels from other cell lines as input variables. The notebook can be found in lb_subs_09_30 folder.

#### Model

We trained xgboost models for each training cell line, and used an average ensemble for the final predictions.

#### Running the code

The preprocessed data was generated using the scripts in create_data_tables folder.  The predictions were made with the scripts in the submissions folder.


#### Installation
Installing and building the programs is straightforward using the Dockerfile provided in the docker folder.

## Discussion 

We have explored the available models to predict TF binding from the sequence, and came to the conclusion, that sometime this one works better sometimes another.
To always use the best inputs possible we incorporated 3 different predictors into the feature table, the baseline, hocomoco and deepbind.
The xgboost trained on the feature table is able to select the best inputs for each TF.


The models had significantly different results when trained, and tested on different cell lines. To make our model more robust, we trained a separate model for each cell line, and used their average value for the final predictions.


## Authors Statement

Each author took part in the discussions where we designed our method. Janos Szalai-Gindl and Sandor Laki  also created the research enviroment, installed and studied the softwares and programs used. Janos-Gindl Szalai studied and reproduced the baseline method provided by the challenge organizers. Akos Rudas explored the dataset. Dezso Ribli created the preprocessing scripts and the submission notebooks. Attila Kiss and Istvan Csabai are the heads of the group and they supervise the work of the other members.

## Acknowledgements

We would like to thank the authors of the baseline method which was provided by the challenge oragizers.

Research was supported by the Novo Nordisk Foundation Interdisciplinary Synergy Programme Grant no.  NNF15OC0016584 and Modernisation and Improvement of Technical Infrastructure for Research and Development Grant no. ITMS 26210120042 .

