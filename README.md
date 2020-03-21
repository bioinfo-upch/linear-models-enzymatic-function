# File: position_analysis.Rmd


The file fluctuation_profiles.R takes four inputs:

- kinetic_parameters.csv: A CSV file containing the relative experimental means of the kinetic parameters.

- electrostatic_potentials.csv: A CSV file containing a dataset with the physicochemical descriptors.

- fluctuations.csv: A CSV file containing a dataset with the stability descriptors derived by normal mode analysis.

- rmsf_backbone.csv: A CSV file containing a dataset with the stability descriptors derived by molecular dynamics.


The outputs are the following figures:

- Comparison of fluctuations, RMSFs, and electrostatic potentials profiles for all the structures under study, highlighting the WT profiles in red.

- Three boxplots of fluctuations, RMSFs, and electrostatic potentials per position along PZAse for structures with a relative k<sub>cat</sub> greater or equal than 50, and other three similar boxplots for structures with a relative k<sub>cat</sub> less than 50.

And CSV files:

- Three CSVs files with the outputs of non-parametric tests per position for the mean value of fluctuations, RMSFs, and electrostatic potentials.

# File: distance_analysis.Rmd


The file distance_analysis takes twelve inputs:

- Activity_models.RData: A set of R objects for the log-linear models for activity

- Efficiency_models.RData: A set of R objects for the log-linear models for efficiency

- Km_models.RData: A set of R objects for the log-linear models for K<sub>M</sub>

- Kcat_models.RData: A set of R objects for the log-linear models for k<sub>cat</sub>

- distance_to_aa_8.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp8

- distance_to_aa_49.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp49

- distance_to_aa_51.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His51

- distance_to_aa_57.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His57

- distance_to_aa_71.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His71

- distance_to_aa_96.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Lys96 

- distance_to_aa_138.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Cys138

- residues.txt: A text file containing the identity of all the residues in PZAse

The output is one figure displaying the 20 nearest residues to each amino-acid of the active site and metal coordination site of wild-type PZAse. Residues associated with selected descriptors are highlighted in blue.
