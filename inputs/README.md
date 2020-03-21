# Input files

## Data

- data_md.csv: A CSV file containing a dataset the kinetic parameters, and the stability (MD-RMSFs) , physicochemical, and geometrical descriptors.
- data_nma.csv: A CSV file containing a dataset the kinetic parameters, and the stability (NMA-fluctuations), physicochemical, and geometrical descriptors.
- descriptors_md.txt: A text file with the real name of all the variables of the dataset with stability descriptors estimated by Molecular Dynamics (MD).
- descriptors_nma.txt: A text file with the real name of all the variables of the dataset with stability descriptors estimated by Normal Mode Analysis (NMA).
- distance_to_aa_8.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp8.
- distance_to_aa_49.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Asp49.
- distance_to_aa_51.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His51.
- distance_to_aa_57.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His57.
- distance_to_aa_71.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and His71.
- distance_to_aa_96.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Lys96.
- distance_to_aa_138.csv: A CSV file with the distances between alpha-carbons of each amino-acid of each structure and Cys138.
- relative_variances.csv: A CSV file containing the relative variances of experimental mesurement for each kinetic parameter
- residues.txt: A text file containing the identity of all the residues in PZAse.

## Structures

This folder contains a total of 36 structures, including the wild-type PZAse crystal 3PL1, and 35 SWISS-MODEL modelled structures of punctual mutations of PZAse for which experimental data is available.

## Models

This folder contains the log-linear models for kcat, Km, efficiency and activity. Each file contains six R objects corresponding to the four individual models for stability (NMA), stabiltiy (MD), physicochemical, and geometrical descriptors, and the two weighted models that differ in the nature of the stability model used in their construction.
