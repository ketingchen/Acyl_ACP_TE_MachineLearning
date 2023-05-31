# Acyl_ACP_TE_MachineLearning
R scripts in this repo are implemented in R version 3.6.
The input files include "TE fatty acid profiles.csv", "TE substrate specificities.csv", and "pairwise comparison of TE sequences encoded with labels.csv".

### The R scripts are used in the order described below for different purposes.
### Step 1: Clustering analysis for fatty acid profiles of each TE --> use script fatty_acid_profile_clustering.R
*Input*: TE fatty acid profiles.csv <br>
*Output*: FA_clustering.rds <br>
fatty_acid_profile_clustering.R depends on two additional scritps: PlotDimRdc.R and PLS.R.
### Step 2: Evaluation of feature importance by random forest (RF) classifier --> use script RF_feature importance.R
*Input*: pairwise comparison of TE sequences encoded with labels.csv; FA_clustering.rds <br>
*Output*: pvalues_10Runs_RF.txt; importance_score_10Runs_RF.txt; importance_rank.txt <br>
<br>
1. Define the instance for RF classifier, response --> comparison of fatty acid profile cluster membership between two TEs (1, same cluster; 0, different clusters); feature --> sequence variation of two TEs at each amino acid position (0, same amino acids; 1, different amino acids)<br>
Step 2.2: construction of RF classifier. The classifier was implemented 10 times using the same dataset to account for the randomness involved in classifier construction.<br>
Step 2.3: calculate the feature importance score and the associted p-values according to the 10 RF classifiers<br>
Step 2.4: calculate the rank of feature importance from high to low <br>
<br>
### Step 3: Further trimming of the results in step 2 to obtain the short list of important features --> use script RF_IFS.R
Incremental feature selection approach was used to identify the RF classifier with the minimum number of features but having an optimal predictive performance. <br>
*Input*: pairwise comparison of TE sequences encoded with labels.csv; TE substrate specificities.csv; importance_rank.txt <br>
*Output*: IFS_evaluation.rds<br>
<br>
### Step 4: Summarize and plot the IFS results --> use script IFS_summary.R
Comparison of MCC was made between every pair of neighbouring models in IFS to determine the minimum number of features needed to reach MCC plateau.<br>
IFS_summary.R requires asterisk.R and errbar.R.<br>
*Input*: IFS_evaluation.rds<br>
*Output*: IFS_model_comparison_mcc.csv<br>
