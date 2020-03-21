########################
## Charging libraries ##
########################

rm(list=ls())
packages <- c('dplyr', 'caret', 'olsrr', 'car', 'rgl','reshape2','ggplot2','biobroom','tidyr',
              'dendextend','matrixStats','cluster','gridExtra','grid','lattice')
for(i in packages){library(package=i,character.only = T)}
source('fpzase.R')

####################################################
## Charging the datasets and statistical modeling ##
####################################################

#Data frame with the mean of kinetic parameters
data.k <- read.table('kinetic_parameters.csv', sep=';', header = T, row.names = 1 )

#Data frame with the variance of kinetic parameters
variances <- read.table('variances.csv', sep=';', header = T, row.names = 1 )
  
#Calculating the mean and variance of the log transformation
data.k = log10((data.k^2)/sqrt(variances+data.k^2))
variances = sqrt(log10((variances/data.k^2)+1))

#Data frame 1: all the physicochemical parameters
data.pc <- read.table('electrostatic_potentials.csv', sep=';', header = T, row.names = 1 )
colnames(data.pc) <- paste('x',1:185,sep='')

#Data frame 1: all the stability NMA parameters
data.nm <- read.table('fluctuations.csv', sep=',', header = T, row.names = 1 )
colnames(data.nm) <- paste('x',186:370,sep='')

#Data frame 1: all the stability MD parameters
data.md <- read.table('rmsf_backbone.csv', sep=';', header = T, row.names = 1 )
colnames(data.md) <- paste('x',4811:4995,sep='')

#Data frame 1: all the geometrical parameters
data.gm <- read.table('geometrical_descriptors.csv', sep=';', header = T, row.names = 1 )
colnames(data.gm) <- paste('x',371:4810,sep='')

#Vector containing the names of all variables
variables <- readLines('descriptors.txt')

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
  data.md.k <- filter.pzase(data.frame(kinetic_param,data.md[mutants,]))
  data.pc.k <- filter.pzase(data.frame(kinetic_param,data.pc[mutants,]))
  data.gm.k <- filter.pzase(data.frame(kinetic_param,data.gm[mutants,]))
  
  #Data frame 3: DF2 with uncorrelated parameters sorted by r2 from simple linear regressions
  # with the kinetic parameter studied 
  data.nm.cor <- correlation.dataset(data.nm.k)$data.cov.red
  data.md.cor <- correlation.dataset(data.md.k)$data.cov.red
  data.pc.cor <- correlation.dataset(data.pc.k)$data.cov.red
  data.gm.cor <- correlation.dataset(data.gm.k)$data.cov.red
  
  ################################
  ## Evaluate regression models ##
  ################################
  dataset = data.nm.cor
  model.nm <- model.construction(data.nm.cor, 10, 6)
  dataset = data.md.cor
  model.md <- model.construction(data.md.cor, 10, 6)
  dataset = data.pc.cor
  model.pc <- model.construction(data.pc.cor, 10, 6)
  dataset = data.gm.cor
  model.gm <- model.construction(data.gm.cor, 10, 6)
  
  ######################
  ## Cross validation ##
  ######################
  
  #Cross validation for selected models. Create 6 folds and for all the combinations,
  #train the model with 5 of them and calculate the RMSE of the remaining fold. Repeat
  #the procedure 1000 times.
  cv.model.nm <- cross.validation(data.nm.k,model.nm,1000,6)
  cv.model.md <- cross.validation(data.md.k,model.md,1000,6)
  cv.model.pc <- cross.validation(data.pc.k,model.pc,1000,6)
  cv.model.gm <- cross.validation(data.gm.k,model.gm,1000,6)
  
  #Create random models, of 6 variables
  rmodel.nm <- random.linear.model(data.nm.k,6)
  rmodel.md <- random.linear.model(data.md.k,6)
  rmodel.pc <- random.linear.model(data.pc.k,6)
  rmodel.gm <- random.linear.model(data.gm.k,6)
  
  #Cross validation for random models. Create 6 folds and for all the combinations,
  #train the model with 5 of them and calculate the RMSE of the remaining fold. Repeat
  #the procedure 1000 times.
  acv.model.nm <- cross.validation(data.nm.k,rmodel.nm,1000,6)
  acv.model.md <- cross.validation(data.md.k,rmodel.md,1000,6)
  acv.model.pc <- cross.validation(data.pc.k,rmodel.pc,1000,6)
  acv.model.gm <- cross.validation(data.gm.k,rmodel.gm,1000,6)
  
  ############################
  ## Model plots and tables ##
  ############################
  
  #Create a correlation heatmap between the variables of the models
  nm.hm <- plot.model.variables(data.nm.cor,names(model.nm$coefficients)[-1],variables)
  md.hm <- plot.model.variables(data.md.cor,names(model.md$coefficients)[-1],variables)
  pc.hm <- plot.model.variables(data.pc.cor,names(model.pc$coefficients)[-1],variables)
  gm.hm <- plot.model.variables(data.gm.cor,names(model.gm$coefficients)[-1],variables)
  
  #Confidence intervals plots for the coefficients of the models
  nm.ci <- gg.intervals(model.nm,variables)
  md.ci <- gg.intervals(model.md,variables)
  pc.ci <- gg.intervals(model.pc,variables,intercept = F)
  gm.ci <- gg.intervals(model.gm,variables,intercept = F)
  
  #Scatter plot for the fitted values of each model and the experimental values
  nm.sp <- plot.scatter.model(model.nm,labels[k])
  md.sp <- plot.scatter.model(model.md,labels[k])
  pc.sp <- plot.scatter.model(model.pc,labels[k])
  gm.sp <- plot.scatter.model(model.gm,labels[k])
  
  #Density plots  of the distribution of RMSE obtained from the cross validation of
  #each model
  cvhist.nm <- cv.rmse.density(cv.model.nm$rmse,acv.model.nm$rmse,
                            model.name = c('Stability Model (NMA)','Random Model'),
                            variable = labels[k])
  cvhist.md <- cv.rmse.density(cv.model.md$rmse,acv.model.md$rmse,
                               model.name = c('Stability Model (MD)','Random Model'),
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
  rownames(model.nm.s) <- c('Stability Model (NMA)','Random Model')
  colnames(model.nm.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','<RMSE>')
  
  model.md <- parameter.name(model.md,variables)
  omodel.md.s <- data.frame(r.squared=round(glance(model.md)$r.squared,3),
                            adj.r.squared=round(glance(model.md)$adj.r.squared,3),
                            p.value=glance(model.md)$p.value,
                            rmse = cv.model.md$mean_rmse)
  amodel.md.s <- data.frame(r.squared=round(glance(rmodel.md)$r.squared,3),
                            adj.r.squared=round(glance(rmodel.md)$adj.r.squared,3),
                            p.value=glance(rmodel.md)$p.value,
                            rmse=acv.model.md$mean_rmse)
  model.md.s <- rbind(omodel.md.s,amodel.md.s)
  model.md.s[,c(1:2,4)] <- round(model.md.s[,c(1:2,4)],4)
  model.md.s[,3] <- formatC(x = model.md.s[,3],format='e',digits = 2)
  rownames(model.md.s) <- c('Stability Model (MD)','Random Model')
  colnames(model.md.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','<RMSE>')

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
  rownames(model.pc.s) <- c('Physicochemical Model','Random Model')
  colnames(model.pc.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','<RMSE>')
  
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
  rownames(model.gm.s) <- c('Geometrical Model','Random Model')
  colnames(model.gm.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','<RMSE>')
  
  model.nm.table <- summary(model.nm)[[4]]
  model.nm.table <- as.data.frame(model.nm.table)
  model.nm.table[,1:3] <- round(model.nm.table[,1:3],4)
  model.nm.table[,4] <- formatC(model.nm.table[,4],format = 'e',digits = 2)
  colnames(model.nm.table) <- c('Coefficient','Std. Error','t-statistic','p-value')
  tbl.nm1 <- tableGrob(model.nm.table, theme = ttheme_default(base_size = 10))
  tbl.nm2 <- tableGrob(model.nm.s,theme = ttheme_default(base_size = 10,
                                                         colhead=list(fg_params = list(parse=T))))
  model.md.table <- summary(model.md)[[4]]
  model.md.table <- as.data.frame(model.md.table)
  model.md.table[,1:3] <- round(model.md.table[,1:3],4)
  model.md.table[,4] <- formatC(model.md.table[,4],format = 'e',digits = 2)
  colnames(model.md.table) <- c('Coefficient','Std. Error','t-statistic','p-value')
  tbl.md1 <- tableGrob(model.md.table, theme = ttheme_default(base_size = 10))
  tbl.md2 <- tableGrob(model.md.s,theme = ttheme_default(base_size = 10,
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
               layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.('Stability Model (NMA)')),gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_stability_nma.tiff',sep = ' '),grob,width=16,height=8,dpi = 300)
  
  grob <- arrangeGrob(tbl.md1, md.hm, tbl.md2, md.sp, md.ci, cvhist.md,
                layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.('Stability Model (MD)')),gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_stability_md.tiff',sep = ' '),grob,width=16,height=8,dpi = 300)
  
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
  w_nma <- data.frame(kinetic_param, Stability=model.nm$fitted.values, Physicochemical=model.pc$fitted.values,
                  Geometrical=model.gm$fitted.values)
  model.wm_nma <- lm(as.formula(paste(colnames(kinetic_param),'~Stability+Physicochemical+Geometrical',sep='')), data=w_nma)
  
  w_md <- data.frame(kinetic_param, Stability=model.md$fitted.values, Physicochemical=model.pc$fitted.values,
                      Geometrical=model.gm$fitted.values)
  model.wm_md <- lm(as.formula(paste(colnames(kinetic_param),'~Stability+Physicochemical+Geometrical',sep='')), data=w_md)
  
  #Scatter plot
  wm_nma.sp <- plot.scatter.model(model.wm_nma,labels[k])
  wm_md.sp <- plot.scatter.model(model.wm_md,labels[k])
  
  #Cross validation for selected models
  cv.model.wm_nma <- cross.validation(w_nma,model.wm_nma,1000,6)
  cvhist.wm_nma <- cv.rmse.density(cv.model.wm_nma$rmse,cv.model.nm$rmse,cv.model.pc$rmse,
                               cv.model.gm$rmse,
                               model.name = c('Weighted Model (NMA)','Stability Model (NMA)',
                                              'Physicochemical Model','Geometrical Model'),
                               variable = labels[k])
  
  cv.model.wm_md <- cross.validation(w_md,model.wm_md,1000,6)
  cvhist.wm_md <- cv.rmse.density(cv.model.wm_md$rmse,cv.model.nm$rmse,cv.model.pc$rmse,
                               cv.model.gm$rmse,
                               model.name = c('Weighted Model (MD)','Stability Model (MD)',
                                              'Physicochemical Model','Geometrical Model'),
                               variable = labels[k])
  
  #Descriptive tables
  model.wm_nma.s <- data.frame(r.squared=round(glance(model.wm_nma)$r.squared,3),
                           adj.r.squared=round(glance(model.wm_nma)$adj.r.squared,3),
                           p.value=glance(model.wm_nma)$p.value,
                           rmse=cv.model.wm_nma$mean_rmse)
  model.wm_nma.s[,c(1:2,4)] <- round(model.wm_nma.s[,c(1:2,4)],4)
  model.wm_nma.s[,3] <- formatC(x = model.wm_nma.s[,3],format='e',digits = 2)
  colnames(model.wm_nma.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','<RMSE>')
  model.wm_nma.s <- rbind(model.wm_nma.s,model.nm.s[1,],model.pc.s[1,],model.gm.s[1,])
  rownames(model.wm_nma.s) <- c('Weighted Model (NMA)','Stability Model (NMA)','Physicochemical Model',
                            'Geometrical Model')
  
  model.wm_nma.table <- summary(model.wm_nma)[[4]]
  model.wm_nma.table <- as.data.frame(model.wm_nma.table)
  model.wm_nma.table[,1:3] <- round(model.wm_nma.table[,1:3],4)
  model.wm_nma.table[,4] <- formatC(model.wm_nma.table[,4],format = 'e',digits = 2)
  colnames(model.wm_nma.table) <- c('Coefficient','Std. Error','t-statistic','p-value')
  
  tbl.wm_nma1 <- tableGrob(model.wm_nma.table, theme = ttheme_default(base_size = 10))
  tbl.wm_nma2 <- tableGrob(model.wm_nma.s, theme = ttheme_default(base_size = 10,
                                                          colhead=list(fg_params = list(parse=T))))
  
  model.wm_md.s <- data.frame(r.squared=round(glance(model.wm_md)$r.squared,3),
                           adj.r.squared=round(glance(model.wm_md)$adj.r.squared,3),
                           p.value=glance(model.wm_md)$p.value,
                           rmse=cv.model.wm_md$mean_rmse)
  model.wm_md.s[,c(1:2,4)] <- round(model.wm_md.s[,c(1:2,4)],4)
  model.wm_md.s[,3] <- formatC(x = model.wm_md.s[,3],format='e',digits = 2)
  colnames(model.wm_md.s) <- c(expression(R^2), bquote(.('Adjusted')~.(expression(R^2)[[1]])),
                            'p-value','<RMSE>')
  model.wm_md.s <- rbind(model.wm_md.s,model.md.s[1,],model.pc.s[1,],model.gm.s[1,])
  rownames(model.wm_md.s) <- c('Weighted Model (MD)','Stability Model (MD)','Physicochemical Model',
                            'Geometrical Model')

  model.wm_md.table <- summary(model.wm_md)[[4]]
  model.wm_md.table <- as.data.frame(model.wm_md.table)
  model.wm_md.table[,1:3] <- round(model.wm_md.table[,1:3],4)
  model.wm_md.table[,4] <- formatC(model.wm_md.table[,4],format = 'e',digits = 2)
  colnames(model.wm_md.table) <- c('Coefficient','Std. Error','t-statistic','p-value')

  tbl.wm_md1 <- tableGrob(model.wm_md.table, theme = ttheme_default(base_size = 10))
  tbl.wm_md2 <- tableGrob(model.wm_md.s, theme = ttheme_default(base_size = 10,
                                                          colhead=list(fg_params = list(parse=T))))
  
  #Put all plots together and save
  lay <- rbind(c(1,1,2,2,2),
               c(1,1,2,2,2),
               c(3,3,2,2,2),
               c(3,3,2,2,2),
               c(4,4,4,4,4),
               c(4,4,4,4,4),
               c(4,4,4,4,4))
  
  grob <- arrangeGrob(tbl.wm_nma1, wm_nma.sp, tbl.wm_nma2, cvhist.wm_nma,
                      layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.(" Weighted Model (NMA)")), gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_weighted_nma.tiff',sep = ' '),grob,width=12,height=8,dpi=300)
  
  grob <- arrangeGrob(tbl.wm_md1, wm_md.sp, tbl.wm_md2, cvhist.wm_md,
               layout_matrix = lay, top = textGrob(bquote(.(labels[k][[1]])~.(" Weighted Model (MD)")), gp=gpar(fontsize=16,font=3)))
  ggsave(paste(colnames(kinetic_param),'_weighted_md.tiff',sep = ' '),grob,width=12,height=8,dpi=300)

  save(model.nm, model.md, model.pc, model.gm, model.wm_nma, model.wm_md,
       file=paste(colnames(data.k[k]),'_models.RData',sep=''))
}