########################
## Charging libraries ##
########################

packages <- c('dplyr', 'caret', 'olsrr', 'car', 'rgl','reshape2','ggplot2','biobroom','tidyr',
              'dendextend','matrixStats','cluster','gridExtra','grid','lattice')
for(i in packages){library(package=i,character.only = T)}
source('fpzase.R')

####################################################
## Charging the datasets and statistical modeling ##
####################################################

#Complete dataset containing all parameters
data <- read.table('finaldataset.csv', sep=';', header = T, row.names = 1 )
#data[1:4] = data[1:4]/100

#Vector containing the names of all variables
variables <- readLines('variables.txt')

#Data frame with the kinetic parameters
data.k <- data[1:4]

#Data frame 1: all the physicochemical parameters
data.pc <- data[5:190]

#Data frame 1: all the stability parameters
data.nm <- data[191:375]

#Data frame 1: all the geometrical parameters
data.gm <- data[376:length(data)]

#Kinetic parameters names to display 
labels <- c(expression(k[cat]),expression(K[m]),expression(Efficiency),expression(Activity))


for(k in 1:length(data.k)){
  
  #########################################################################
  ## Selecting the most significant variables of each variable's cluster ##
  #########################################################################
  
  #Select the mutants that have all the kinetic parameters measured 
  mutants <- data.k[,k] != 0 & !is.na(data.k[,k])
  
  #Create a data frame with just one of the kinetic parameters
  kinetic_param <- data.frame(data.k[mutants,k])
  colnames(kinetic_param) <- colnames(data.k[k])
  
  #Data frame 2: DF1 + kinetic parameter of studied - mutants with no data for that kinetic
  # - parameters with missing data or variance equal to 0
  data.nm.k <- filter.pzase(data.frame(kinetic_param,data.nm[mutants,]))
  data.pc.k <- filter.pzase(data.frame(kinetic_param,data.pc[mutants,]))
  data.gm.k <- filter.pzase(data.frame(kinetic_param,data.gm[mutants,]))
  
  #Data frame 3: DF2 with uncorrelated parameters sorted by r2 from simple linear regressions
  # with the kinetic parameter studied 
  data.nm.cor <- correlation.dataset(data.nm.k)$data.cov.red
  data.pc.cor <- correlation.dataset(data.pc.k)$data.cov.red
  data.gm.cor <- correlation.dataset(data.gm.k)$data.cov.red
  
  ################################
  ## Evaluate regression models ##
  ################################
  
  model.nm <- model.construction(data.nm.cor, 10, 6)
  model.pc <- model.construction(data.pc.cor, 10, 6)
  model.gm <- model.construction(data.gm.cor, 10, 6)
  
  ######################
  ## Cross validation ##
  ######################
  
  #Cross validation for selected models. Create 6 folds and for all the combinations,
  #train the model with 5 of them and calculate the RMSE of the remaining fold. Repeat
  #the procedure 1000 times.
  cv.model.nm <- cross.validation(data.nm.k,model.nm,1000,6)
  cv.model.pc <- cross.validation(data.pc.k,model.pc,1000,6)
  cv.model.gm <- cross.validation(data.gm.k,model.gm,1000,6)
  
  #Create random models, of 6 variables
  rmodel.nm <- random.linear.model(data.nm.k,6)
  rmodel.pc <- random.linear.model(data.pc.k,6)
  rmodel.gm <- random.linear.model(data.gm.k,6)
  
  #Cross validation for random models. Create 6 folds and for all the combinations,
  #train the model with 5 of them and calculate the RMSE of the remaining fold. Repeat
  #the procedure 1000 times.
  acv.model.nm <- cross.validation(data.nm.k,rmodel.nm,1000,6)
  acv.model.pc <- cross.validation(data.pc.k,rmodel.pc,1000,6)
  acv.model.gm <- cross.validation(data.gm.k,rmodel.gm,1000,6)
  
  ############################
  ## Model plots and tables ##
  ############################
  
  #Create a correlation heatmap between the variables of the models
  nm.hm <- plot.model.variables(data.nm.cor,names(model.nm$coefficients)[-1],variables)
  pc.hm <- plot.model.variables(data.pc.cor,names(model.pc$coefficients)[-1],variables)
  gm.hm <- plot.model.variables(data.gm.cor,names(model.gm$coefficients)[-1],variables)
  
  #Confidence intervals plots for the coefficients of the models
  nm.ci <- gg.intervals(model.nm,variables)
  pc.ci <- gg.intervals(model.pc,variables,intercept = F)
  gm.ci <- gg.intervals(model.gm,variables,intercept = F)
  
  #Scatter plot for the fitted values of each model and the experimental values
  nm.sp <- plot.scatter.model(model.nm,labels[k])
  pc.sp <- plot.scatter.model(model.pc,labels[k])
  gm.sp <- plot.scatter.model(model.gm,labels[k])
  
  #Density plots  of the distribution of RMSE obtained from the cross validation of
  #each model
  cvhist.nm <- cv.rmse.density(cv.model.nm$rmse,acv.model.nm$rmse,
                            model.name = c('Stability Model','Random Model'),
                            variable = labels[k])
  cvhist.pc <- cv.rmse.density(cv.model.pc$rmse,acv.model.pc$rmse,
                            model.name = c('Physicochemical Model','Random Model'),
                            variable = labels[k])
  cvhist.gm <- cv.rmse.density(cv.model.gm$rmse,acv.model.gm$rmse,
                            model.name = c('Geometrical Model','Random Model'),
                            variable = labels[k])
  
  #Descriptive tables
  model.nm <- parameter.name(model.nm,variables)
  omodel.nm.s <- data.frame(r.squared=round(glance(model.nm)$r.squared,3),
                           adj.r.squared=round(glance(model.nm)$adj.r.squared,3),
                           p.value=glance(model.nm)$p.value,
                           rmse = cv.model.nm$mean_rmse)
  amodel.nm.s <- data.frame(r.squared=round(glance(rmodel.nm)$r.squared,3),
                            adj.r.squared=round(glance(rmodel.nm)$adj.r.squared,3),
                            p.value=glance(rmodel.nm)$p.value,
                            rmse=acv.model.nm$mean_rmse)
  model.nm.s <- rbind(omodel.nm.s,amodel.nm.s)
  model.nm.s[,c(1:2,4)] <- round(model.nm.s[,c(1:2,4)],4)
  model.nm.s[,3] <- formatC(x = model.nm.s[,3],format='e',digits = 2)
  rownames(model.nm.s) <- c('Stability Model','Random model')
  colnames(model.nm.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','mean RMSE')

  model.pc <- parameter.name(model.pc,variables)
  omodel.pc.s <- data.frame(r.squared=round(glance(model.pc)$r.squared,3),
                           adj.r.squared=round(glance(model.pc)$adj.r.squared,3),
                           p.value=glance(model.pc)$p.value,
                           rmse=cv.model.pc$mean_rmse)
  amodel.pc.s <- data.frame(r.squared=round(glance(rmodel.pc)$r.squared,3),
                            adj.r.squared=round(glance(rmodel.pc)$adj.r.squared,3),
                            p.value=glance(rmodel.pc)$p.value,
                            rmse=acv.model.pc$mean_rmse)
  model.pc.s <- rbind(omodel.pc.s,amodel.pc.s)
  model.pc.s[,c(1:2,4)] <- round(model.pc.s[,c(1:2,4)],4)
  model.pc.s[,3] <- formatC(x = model.pc.s[,3],format='e',digits = 2)
  rownames(model.pc.s) <- c('Physicochemical Model','Random model')
  colnames(model.pc.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','mean RMSE')
  
  model.gm <- parameter.name(model.gm,variables)
  omodel.gm.s <- data.frame(r.squared=round(glance(model.gm)$r.squared,3),
                           adj.r.squared=round(glance(model.gm)$adj.r.squared,3),
                           p.value=glance(model.gm)$p.value,
                           rmse=cv.model.gm$mean_rmse)
  amodel.gm.s <- data.frame(r.squared=round(glance(rmodel.gm)$r.squared,3),
                            adj.r.squared=round(glance(rmodel.gm)$adj.r.squared,3),
                            p.value=glance(rmodel.gm)$p.value,
                            rmse=acv.model.gm$mean_rmse)
  model.gm.s <- rbind(omodel.gm.s,amodel.gm.s)
  model.gm.s[,c(1:2,4)] <- round(model.gm.s[,c(1:2,4)],4)
  model.gm.s[,3] <- formatC(x = model.gm.s[,3],format='e',digits = 2)
  rownames(model.gm.s) <- c('Geometrical Model','Random model')
  colnames(model.gm.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','mean RMSE')
  
  model.nm.table <- summary(model.nm)[[4]]
  model.nm.table <- as.data.frame(model.nm.table)
  model.nm.table[,1:3] <- round(model.nm.table[,1:3],4)
  model.nm.table[,4] <- formatC(model.nm.table[,4],format = 'e',digits = 2)
  colnames(model.nm.table) <- c('Coefficient','Std. Error','t-statistic','p-value')
  tbl.nm1 <- tableGrob(model.nm.table, theme = ttheme_default(base_size = 10))
  tbl.nm2 <- tableGrob(model.nm.s,theme = ttheme_default(base_size = 10,
                                                         colhead=list(fg_params = list(parse=T))))
  
  model.pc.table <- summary(model.pc)[[4]]
  model.pc.table <- as.data.frame(model.pc.table)
  model.pc.table[,1:3] <- round(model.pc.table[,1:3],4)
  model.pc.table[,4] <- formatC(model.pc.table[,4],format = 'e',digits = 2)  
  colnames(model.pc.table) <- c('Coefficient','Std. Error','t-statistic','p-value')
  tbl.pc1 <- tableGrob(model.pc.table, theme = ttheme_default(base_size = 10))
  tbl.pc2 <- tableGrob(model.pc.s,theme = ttheme_default(base_size = 10,
                                                         colhead=list(fg_params = list(parse=T))))
  model.gm.table <- summary(model.gm)[[4]]
  model.gm.table <- as.data.frame(model.gm.table)
  model.gm.table[,1:3] <- round(model.gm.table[,1:3],4)
  model.gm.table[,4] <- formatC(model.gm.table[,4],format = 'e',digits = 2)
  colnames(model.gm.table) <- c('Coefficient','Std. Error','t-statistic','p-value')
  tbl.gm1 <- tableGrob(model.gm.table, theme = ttheme_default(base_size = 10))
  tbl.gm2 <- tableGrob(model.gm.s,theme = ttheme_default(base_size = 10,
                                                         colhead=list(fg_params = list(parse=T))))
  #Put all plots together and save
  lay <- rbind(c(1,1,1,4,4,4,4,2,2,2),
               c(1,1,1,4,4,4,4,2,2,2),
               c(3,3,3,4,4,4,4,2,2,2),
               c(5,5,5,5,5,6,6,6,6,6),
               c(5,5,5,5,5,6,6,6,6,6),
               c(5,5,5,5,5,6,6,6,6,6))
  
  grob <- arrangeGrob(tbl.nm1, nm.hm, tbl.nm2, nm.sp, nm.ci, cvhist.nm,
               layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.('Stability Model')),gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_stability.tiff',sep = ' '),grob,width=16,height=8,dpi = 300)
  
  grob <- arrangeGrob(tbl.pc1, pc.hm, tbl.pc2, pc.sp, pc.ci, cvhist.pc,
               layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.('Physicochemical Model')),gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_physicochemical.tiff',sep = ' '),grob,width=16,height=8,dpi=300)
  
  grob <- arrangeGrob(tbl.gm1, gm.hm, tbl.gm2, gm.sp, gm.ci, cvhist.gm,
               layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.('Geometrical Model')),gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_geometric.tiff',sep = ' '),grob,width=16,height=8,dpi=300)

  #############################
  ## Create a weighted model ##
  #############################
  
  #Defining the model
  w <- data.frame(kinetic_param, Stability=model.nm$fitted.values, Physicochemical=model.pc$fitted.values,
                  Geometrical=model.gm$fitted.values)
  model.wm <- lm(as.formula(paste(colnames(kinetic_param),'~Stability+Physicochemical+Geometrical',sep='')), data=w)
  
  #Scatter plot
  wm.sp <- plot.scatter.model(model.wm,labels[k])
  
  #Cross validation for selected models
  cv.model.wm <- cross.validation(w,model.wm,1000,6)
  cvhist.wm <- cv.rmse.density(cv.model.wm$rmse,cv.model.nm$rmse,cv.model.pc$rmse,
                               cv.model.gm$rmse,
                               model.name = c('Weighted Model','Stability Model',
                                              'Physicochemical Model','Geometrical Model'),
                               variable = labels[k])
  
  #Descriptive tables
  model.wm.s <- data.frame(r.squared=round(glance(model.wm)$r.squared,3),
                           adj.r.squared=round(glance(model.wm)$adj.r.squared,3),
                           p.value=glance(model.wm)$p.value,
                           rmse=cv.model.wm$mean_rmse)
  model.wm.s[,c(1:2,4)] <- round(model.wm.s[,c(1:2,4)],4)
  model.wm.s[,3] <- formatC(x = model.wm.s[,3],format='e',digits = 2)
  colnames(model.wm.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','mean RMSE')
  model.wm.s <- rbind(model.wm.s,model.nm.s[1,],model.pc.s[1,],model.gm.s[1,])
  rownames(model.wm.s) <- c('Weighted Model','Stability model','Physicochemical model',
                            'Geometrical model')


  
  model.wm.table <- summary(model.wm)[[4]]
  model.wm.table <- as.data.frame(model.wm.table)
  model.wm.table[,1:3] <- round(model.wm.table[,1:3],4)
  model.wm.table[,4] <- formatC(model.wm.table[,4],format = 'e',digits = 2)
  colnames(model.wm.table) <- c('Coefficient','Std. Error','t-statistic','p-value')

  tbl.wm1 <- tableGrob(model.wm.table, theme = ttheme_default(base_size = 10))
  tbl.wm2 <- tableGrob(model.wm.s, theme = ttheme_default(base_size = 10,
                                                          colhead=list(fg_params = list(parse=T))))
  
  #Put all plots together and save
  lay <- rbind(c(1,1,2,2,2),
               c(1,1,2,2,2),
               c(3,3,2,2,2),
               c(3,3,2,2,2),
               c(4,4,4,4,4),
               c(4,4,4,4,4),
               c(4,4,4,4,4))
  grob <- arrangeGrob(tbl.wm1, wm.sp, tbl.wm2, cvhist.wm,
               layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.(" Weighted Model")), gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_weighted.tiff',sep = ' '),grob,width=12,height=8,dpi=300)

  save(model.nm,model.pc,model.gm,model.wm,
       file=paste(colnames(data.k[k]),'_models.RData',sep=''))
}