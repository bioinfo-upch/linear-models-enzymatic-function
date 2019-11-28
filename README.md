

# File: fpzase.R 


The file fpzase.R contain all the functions for the statistical analysis


# File: fluctuation_profiles.R


The file fluctuation_profiles.R takes one input:

- finaldataset.csv: A CSV file containing a dataset with the kinetic parameters, and the stability, physicochemical, and geometrical descriptors

The outputs are the following figures:

- Comparison of fluctuation profiles for two mutated and the wild-type PZAse

- A boxplot of fluctuation per position for the flap region of PZAse (His51 - His71). Each position has two boxplots, one for PZAses with catalytic constant under the median value and the other above it

- A boxplot of fluctuation per pesition along all the PZAse for enzymes with a catalytic constant above the median value

- A boxplot of fluctuation per pesition along all the PZAse for enzymes with a catalytic constant under the median value


# File: models_construction.R


The file models_construction.R takes two inputs:

- finaldataset.csv: A CSV file containing a dataset with the kinetic parameters, and the stability, physicochemical, and geometrical descriptors 

- variables.txt: A text file with the real name of all the variables of the dataset

The outputs are a group of figures for each modeled kinetic parameter. Three of them are for each individual model (stability, physicochemical and geometrical model) that contains:

- A table with a summary of the model

- A table comparing the prediction RMSE and R2 of the model against a random model of the same kind of predictors

- A scatter plot between fitted and experimental values for the model

- A heatmap for the pearson correlation coefficient between the variables of the model

- A plot of the confidence intervals for the coefficients of the variables of the model

- A comparison of the density plots for the RMSE of prediction of the model and a random model of the same kind of predictors

The fourth plot is for a weighted model fitted with the predictions of the individual models that contains:

- A table with a summary of the model

- A table comparing the prediction RMSE and R2 of the weighted model and the individual models

- A scatter plot between fitted and experimental values for the model

- A comparison of the density plots for the RMSE of prediction between the weighted model and the individual models


# File: distance_matrix_analysis.R


The file distance_matrix_analysis takes nine inputs:

- distance_to_aa_8.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp8

- distance_to_aa_49.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp49

- distance_to_aa_51.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His51

- distance_to_aa_57.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His57

- distance_to_aa_71.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His71

- distance_to_aa_96.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Lys96 

- distance_to_aa_138.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Cys138

- residues.txt: A text file containing the identity of all the residues in PZAse

- selected_residues.txt: A txt file with the names of the residues associated with the selected parameters of the models constructed by models_construction.R

The output is one figure displaying the 20 nearest residues to each amino-acid of the active site and metal coordination site of wild-type PZAse. Residues associated with selected descriptors are highlighted in blue.
