# File: position_analysis.Rmd


The file fluctuation_profiles.R takes one input:

- finaldataset.csv: A CSV file containing a dataset with the kinetic parameters, and the stability, physicochemical, and geometrical descriptors

kinetic_parameters.csv: A CSV file containing the relative experimental means of the kinetic parameters.

electrostatic_potentials.csv: A CSV file containing a dataset with the physicochemical descriptors.

fluctuations.csv: A CSV file containing a dataset with the stability descriptors derived by normal mode analysis.

rmsf_backbone.csv: A CSV file containing a dataset with the stability descriptors derived by molecular dynamics.



The outputs are the following figures:

- Comparison of fluctuation profiles for two mutated and the wild-type PZAse

- A boxplot of fluctuation per position for the flap region of PZAse (His51 - His71). Each position has two boxplots, one for PZAses with catalytic constant under the median value and the other above it

- A boxplot of fluctuation per pesition along all the PZAse for enzymes with a catalytic constant above the median value

- A boxplot of fluctuation per pesition along all the PZAse for enzymes with a catalytic constant under the median value


# File: distance_analysis.Rmd


The file distance_matrix_analysis takes nine inputs:

- Activity_models.RData

- Efficiency_models.RData

- Km_models.RData

- Kcat_models.RData

- distance_to_aa_8.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp8

- distance_to_aa_49.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp49

- distance_to_aa_51.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His51

- distance_to_aa_57.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His57

- distance_to_aa_71.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His71

- distance_to_aa_96.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Lys96 

- distance_to_aa_138.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Cys138

- residues.txt: A text file containing the identity of all the residues in PZAse

The output is one figure displaying the 20 nearest residues to each amino-acid of the active site and metal coordination site of wild-type PZAse. Residues associated with selected descriptors are highlighted in blue.
