######################################################################
## Eliminate covariates with NA observations or variance equal to 0 ##
######################################################################

filter.pzase <- function(data){
  # Define an empty vector.
  d <- c()
  # For each column in the dataset, if there is a missing value, all the observations
  # are equal to zero, or the standard deviation is zero, add that variable to the
  # empty vector.
  for(i in 1:ncol(data)){
    if(all(!is.na(data[i])) & !(all(data[i]==0)) & (sd(data[,i])!=0)){
      next
    } else {
      d <- c(d,i)
      next
    }
  }
  # If there is no variables in the empty vector to filter, display it, otherwise 
  # remove them from the dataset.
  if(length(d)==0){
    print('no variables filtered!!!')
    data <- data
    }
  else{data <- data[,-c(d)]}
}

###############################
## Create a reduced data set ##
###############################

correlation.dataset <- function(data,cutoff=0.8){
  
  # Inputs:
  # data = dataset to reduce. It must include the variable to model in the first column
  # cutoff = cutoff for the pearson correlation coefficient between variables
  
  # Outputs:
  # data.cov.red = reduced dataset with variables ordered according to the p-value of 
  #                simple linear regressions and with |R| between them less than the cutoff
  # pval.order = dataset with the p-value of linear regression for all the variables
  # clusters = clusters of correlation (not used)
  
  # Vector with the name of the variables to evaluate. Do not include the first column
  # (variable to model)
  varnames <- colnames(data[-1])
  # Create a dataset without the variable to model and define an empty vector
  data.cov <- data[-1]
  pv <- c()
  
  # Correlate each variable with the variable of interest and store the p-value of the
  # linear regression in the vector created above
  for(var in varnames){
    formula <- as.formula(paste(colnames(data[1]),'~',var,sep=''))
    reg <- lm(formula,data)
    p <- glance(reg)[[5]]
    pv <- c(pv,p)
  }
  # Order the p-values in ascending order
  pvmo <- order(pv)
  
  # Reorder the variables in the dataset according to the regression p-value
  varnames <- varnames[pvmo] 
  data.cov.ord <- data.cov[pvmo]
  
  # Create a vector to store the name of the variables of the reduced dataset
  surv <- varnames
  clusters <- list()
  # For each variable in the order dataset, calculate the aboslute valur of the
  # pearson correlation coefficient of it with the rest of the variables, then 
  # eliminate the variables with a value higher than the cutoff. Repeat the same
  # with the variables that remain in the dataset.
  for(i in 1:length(varnames)){
    if(varnames[i] %in% surv){
      p <- cor(data.cov.ord[i],data.cov.ord[-i])
      p2 <- abs(as.numeric(p))
      n <- colnames(p)
      kill <- (p2 > cutoff)
      kill <- n[kill]
      clusters[[varnames[i]]]= kill
      surv <- setdiff(surv,kill)
    }
  }
  # Store the pvalue of all the variables
  pval.order <- data.frame(variable = pvmo, pval = sort(pv))
  # Create a reduced dataset with variables ordered according to the p-value they
  # produce in simple linear regressions and with a absolute pearson correlation 
  # coefficient between them less than the cutoff. 
  data.cov.red = data.frame(data[1],data.cov[,surv])
  
  return(list(data.cov.red = data.cov.red,pval.order = pval.order,clusters = clusters))
}

#######################################
## Models construction and selection ##
#######################################

# Stepwise Akaike Information Criterion Forward Regression 

staic <- function(data, maxk=1, prev=numeric(0)){
  
  # Inputs
  # data = dataset with the variables to select
  # maxk = number of variables to include in the model
  # prev = empty vector to store the indexes of the selected variables
  
  # Output
  # A character vector with the name of the selected variables
  
  ic <- numeric(0)
  vid<- numeric(0)
  
  # For each variable not included in the vector prev, calculate the Akaike
  # Information Criterion of the model for that variable and the variables in
  # the vector prev.
  for(i in setdiff(1:(ncol(data)-1),prev)){
    formula <- as.formula(paste(colnames(data[1]),'~.',sep=''))
    ml <- lm(formula, data = data[,c(1,1+c(prev,i))])
    ic <- c(ic,AIC(ml))
    vid <- c(vid,i)
  }
  # Add to the model the added variable that produced the minimum AIC. Repeat the 
  # procedure until you get the number of variables specified.
  nxt <- c(prev,vid[which(ic==min(ic))])
  if(length(nxt) >= maxk){
    return(colnames(data[nxt+1]))
  } else{
    staic(data, maxk, nxt)
  }
}

# Construction of linear models

model.construction <- function(data,n,m){
  
  # Inputs
  # data: dataset with the variables to select
  # n = number of variables to select with the Stepwise AIC Forward Regression 
  # m = number of variables to include in the model
  
  # Output
  # A linear model with the selected predictors
  
  # Select variables using the Stepwise AIC Forward Regression
  variables <- staic(data, maxk=n)
  formula <- as.formula(paste(colnames(data[1]),'~',
                              paste(variables,collapse = '+'),sep=''))
  model <- lm(formula,data)
  regs <- ols_step_all_possible(model)
  bestmodel <- c()
  index <- 0
  # Select the model among all the combinations of the previously selected variables
  # with the highest adjuted R-squared and with a number of variables equal to the input m
  for(i in strsplit(regs$predictors, ' ')){
    index <- index+1
    if(length(bestmodel)>=1){
      break
    }
    if(length(i)==m){
      bestmodel <- c(bestmodel,index)
    }
  }
  # Create a model object as the function output
  selected.variables <- paste(strsplit(regs[bestmodel,3],' ')[[1]],collapse='+')
  formula <- as.formula(paste(colnames(data[1]),'~',
                              selected.variables,sep=''))
  final.model <- lm(formula,data)
  
  return(final.model)
}


# Generate a random linear model for a given data set
random.linear.model <- function(data,n){
  # Inputs:
  # data = dataset with a variable to model and a set of independent variables
  # n = number of variables of the random model
  
  # Output:
  # A random linear model 
  variables <- sample(colnames(data[-1]),n)
  formula <- as.formula(paste(colnames(data[1]),'~',
                              paste(variables,collapse = '+'),sep=''))
  model <- lm(formula,data)
  return(model)
}

##################
## Heatmap Plot ##
##################

plot.model.variables <- function(data,variables,realnames){
  # Inputs:
  # data = dataset with variables of interest
  # variables = character vector with the name of the variables to consider
  # realnames = character vector with the real name of the variables 
  
  # Output:
  # A heatmap between the selected variables
  
  # Remove the variable to model and create a correlation matrix between
  # the selected variables
  data_cov.cor <- data[variables]
  index <- as.numeric(unlist(strsplit(paste(variables,collapse = ''),split='x'))[-1])
  colnames(data_cov.cor) <- realnames[index] 
  data_cov.cor <- cor(data_cov.cor)
  melted_cormat <- melt(data_cov.cor)
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
  }
  upper_tri <- get_upper_tri(data_cov.cor)
  
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  # Heatmap
  ggheatmap<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed()

  ggheatmap + 
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
}

##########################
## Confidence intervals ##
##########################

gg.intervals <- function(model,variables,intercept=T){
  # Inputs
  # model = model object 
  # variables = character vector with the real names of all the variables
  # intercept = logical value indicating wheter the intercept coefficient should be plotted
  
  # Output:
  # A plot displaying the value of the coefficient the width of their confident intervals
  
  #Confidence intervals for coefficients  
  td <- tidy(model, conf.int = TRUE)
  colnames(td)[1] <- 'Parameters'
  index <- as.numeric(unlist(strsplit(paste(names(model$coefficients)[-1]
                                            ,collapse = ''),split='x'))[-1])
  if(!intercept){
    td <- td[-1,]
    td[,1] <- variables[index]
  } else {
    td[-1,1] <- variables[index]
  }
  td$Parameters <- factor(td$Parameters, levels=td$Parameters)
  ggplot(td, aes(estimate, Parameters, color=Parameters)) +
    geom_point() + geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) + ylab('Coefficients') +
    xlab('95% Confidence Intervals')
}

########################
## Model scatter plot ##
########################

plot.scatter.model <- function(model,kinetic_label){
  # Inputs:
  # model = model object whose predictions will be plotted against the real values
  # kinetic_label = name of the modeled variable to show in the plot
  
  # Output:
  # A scatter plot between fitted and real values
  
  kinetic_parameter <- names(model$model[1])
  aug <- augment(model, conf.int = TRUE)
  pred <- predict(model)
  aug <- cbind(aug,pred)
  ggplot(aug, aes_string(x=kinetic_parameter,y='.fitted',label='rownames(aug)')) + geom_point() +
    geom_abline(col=2) + ylab(bquote(.("Predicted WT normalized")~.(kinetic_label[[1]]))) +
    xlab(bquote(.("WT normalized")~.(kinetic_label[[1]]))) + geom_text(hjust=-0.25,size=2.5)
}

######################
## Cross validation ##
######################

cross.validation <- function(data,model,iterations,nfolds){
  # Inputs:
  # data = dataset with observations to create the folds
  # model = model object with variables to validate
  # iterations = number of times to perform cross validation
  # nfolds = number of folds in which the dataset will be splitted
  
  # Outputs:
  # coefficients = estimated coefficients for each fitting.
  # mean_rmse = mean root mean square error over all the model validations.
  # rmse = vector containing the root mean square errors for each model evaluation. 
   
  # Create a matrix to store the coefficients estimated in each iteration
  coefficients <- matrix(nrow=0,ncol=length(model$coefficients))
  # Create a vector to store the root mean square error of prediction for each iteration
  rmse <- c()
  # For each iteration, create n folds. Then for all the combinations of folds, fit a 
  # model with n-1 folds and test the model with the not used fold. Store the rmse for
  # the fold used to test.
  for(i in 1:iterations){
    folds <- createFolds(1:nrow(data),nfolds)
    for(j in 1:length(folds)){
      formula <- paste(names(model$coefficients)[-1],collapse = '+')
      formula <- as.formula(paste(names(data[1]),'~',formula, sep=''))
      samplemodel <- lm(formula, data[-folds[[j]],])
      coefficients <- rbind(coefficients,samplemodel$coefficients)
      test <- predict(object = samplemodel, data[folds[[j]],])
      rmse <- c(rmse,sqrt(sum((test - data[folds[[j]],1])^2)/length(test)))
    }
  }
  return(list(coefficients=coefficients, mean_rmse = mean(rmse), rmse=rmse))
}
  
cv.rmse.density <- function(..., model.name, variable){
  # Inputs:
  # ...: Vectors of rmse 
  # model.name = character vector with the name of the models which evaluation produce
  #              the provided rmse vectors
  # variable: expression that indicated the name of the modeled variable
  
  # Output:
  # Density plots that shown the distribution of rmse for different croos validated models
  
  data <- list(...)
  fill <- c()
  RMSE <- c()
  for(i in 1:length(data)){
    fill <- c(fill,rep(model.name[i],length(data[[i]])))
    RMSE <- c(RMSE,data[[i]])
  }
  fill <- factor(x = fill, levels = unique(fill))
  dat = data.frame(RMSE = RMSE, Model = fill)
  plot <- ggplot(dat, aes(x=RMSE, fill=Model)) + geom_density(alpha=0.5) +
    theme(axis.title = element_text(size = 10), axis.text = element_text(20)) +
    labs(x=bquote(.("RMSE distribution of predicted WT normalized ")~.(variable[[1]])),
         y='Density')
  return(plot)
}                


#########################################################
## Rename the selected variables for plotting purposes ##
#########################################################

parameter.name <- function(model,variables){
  # Inputs:
  # model = model object whith variables to rename
  # variables = character vector with the real name of all the variables studied
  
  # Output: 
  # A model object with renamed variables.
  
  index <- unlist(strsplit(paste(rownames(summary(model)[[4]])[-1],collapse=''),'x'))
  index <- index[-1]
  names <- variables[as.numeric(index)]
  names(model$coefficients)[-1] <- names
  return(model)
}
