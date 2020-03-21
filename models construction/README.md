
Author: Rydberg Supo Escalante

# File: fpzase.R 


The file fpzase.R contain all the functions for the statistical analysis

- filter.pzase: function to filter columns with missing values or variance equal to zeto.

- correlation.dataset: function to reduce a dataset by sorting ascendingly its covariables according to the p-value of simple linear regressions against a variable of interest. Variables with a correlation coefficent greater than |0.8| are dropped and only the variable with the smallest p-value is conserved.

- staic: function to perform a stepwise Akaike Information Criterion forward regression in a given dataset and generate the best model for a fixed number of variables

- model.construction: function to actually construct the best predictive model for a given dataset. It calls the staic function to generate a subset of covariables and evaluate all the models for a fixed number of variables between them. The model with the best $R^2$ is selected.

- random.linear.model

- plot.model.variables

- gg.intervals

- plot.scatter.model

- cross.validation

- cv.rmse.density

- parameter.name


# File: models_construction.R


The file models_construction.R takes seven inputs:

- kinetic_parameters.csv: A CSV file containing the relative experimental means of the kinetic parameters.

- relative_variancs.csv: A CSV file containing the relative experimental variances of the kinetic parameters.

- electrostatic_potentials.csv: A CSV file containing a dataset with the physicochemical descriptors.

- fluctuations.csv: A CSV file containing a dataset with the stability descriptors derived by normal mode analysis.

- rmsf_backbone.csv: A CSV file containing a dataset with the stability descriptors derived by molecular dynamics.

- geometrical_descriptors.csv: A CSV file containing a dataset with the geometrical descriptors.

- descriptors.txt: A text file with the real name of all the descriptros for each data set

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
