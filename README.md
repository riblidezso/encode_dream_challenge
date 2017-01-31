# Eotvos Lorand University team, ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge

Authors: Dezso Ribli, Janos Szalai-Gindl, Akos Rudas, Sandor Laki, Attila Kiss, Istvan Csabai

Affiliations: Eotvos Lorand University

The submission can be made public.

## Background/Intro

We joined very late to the challenge, and our submission in based on the methods of the baseline script. We are hoping to develop a more advanced model for the next final submission round.


## Methods

#### Features

We have used the following variables:
- Motif scores as calculated in the baseline method. Used parallelism to speed up things.
- Dnase raw fold coverage values for every cell line in the position. Also parallelized.
- Chipseq label data in other cell lines at the position if available. The regions in these files did not exactly match the ones described in test_regions.blacklistfiltered.bed.gz which is the file that had to be used for the submissions. So we made extended label files which contain these regions, by filling in missing values with zeroes.
- The predicted labels were transformed as the following (U->0,A->0,B->1)

Please note that our leaderboard submission is different from the final submissions, as it does not use the chipseq labels from other cell lines as input variables. The notebook can be found in lb_subs_09_30 folder.

#### Model
We have used logistic regression, trained in a stochastic gradient fashion, using scikit-learn to predict the labels.

#### Running the code
The preprocessed data was generated using the scripts in create_data_tables folder.  The predictions were made in jupyter notebooks and these notebooks can be found in the final_subs_09_30 folder. To run the notebooks just click run all.  Please note that the notebooks and the scripts were developed in parallelel to running them, because of the strict deadline, therefore some very small modification might be needed to run them succesfully.

#### Installation
Installing and building the programs is very straightforward using the Dockerfile provided in the docker folder.



## Discussion 

As stated before we joined very late to the challenge. By using the chipseq label data from other cell lines, our predictions are much better for the regions where train data was given. This makes our submissions weak on the test regions, but stronger in a real world screnario descirbed in the challenge description. We are hoping to develop a more advanced model for the next final submission round. 

## Authors Statement

Each author took part in the discussions where we designed our method. Janos Szalai-Gindl and Sandor Laki  also created the research enviroment, installed and studied the softwares and programs used. Janos-Gindl Szalai studied and reproduced the baseline method provided by the challenge organizers. Akos Rudas explored the dataset. Dezso Ribli created the preprocessing scripts and the submission notebooks. Attila Kiss and Istvan Csabai are the heads of the group and they supervise the work of the other members.

## Acknowledgements

We would like to thank the authors of the baseline method which was provided by the challenge oragizers.

Research was supported by the Novo Nordisk Foundation Interdisciplinary Synergy Programme Grant no.  NNF15OC0016584 and Modernisation and Improvement of Technical Infrastructure for Research and Development Grant no. ITMS 26210120042 .

