# Acyl_ACP_TE_MachineLearning
The R codes used to select important residues impacting TE substrate specificities are provided in "TE_randomForest.R". R codes are implemented in R version 3.6.
The input files include "TE fatty acid profiles.csv" and "pairwise comparison of TE sequences encoded with labels.csv".

###### The following R scripts are used following the order described below for different purposes
## Step 1: Clustering analysis for fatty acid profiles of each TE --> use script fatty_acid_profile_clustering.R
## Step 2: Evaluation of feature importance by random forest (RF) classifier --> use script RF_feature importance.R
Step 2.1: define the instance for RF classifier, response --> comparison of fatty acid profile cluster membership between two TEs (1, same cluster; 0, different clusters); feature --> sequence variation of two TEs at each amino acid position (0, same amino acids; 1, different amino acids)<br>
Step 2.2: construction of RF classifier. The classifier was implemented 10 times using the same dataset to account for the randomness involved in classifier construction.<br>
Step 2.3: calculate the feature importance score and the associted p-values according to the 10 RF classifiers<br>
Step 2.4: calculate the rank of feature importance from high to low <br>
## Step 3: Further trimming of the results in step 2 to obtain the short list of important features. Incremental feature selection approach was used to identify the RF classifier with the minimum number of features but having an optimal predictive performace --> use script 