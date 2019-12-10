########################
## Charging libraries ##
########################

packages <- c('dplyr', 'caret', 'olsrr', 'car', 'rgl','reshape2','ggplot2','biobroom','tidyr',
              'dendextend','matrixStats','cluster','gridExtra','grid','lattice')
for(i in packages){library(package=i,character.only = T)}


##########################
## Loading the datasets ##
##########################

#Complete dataset containing all parameters
data <- read.table('finaldataset.csv', sep=';', header = T, row.names = 1 )

#Data frame with the kinetic parameters
data.k <- data[1:4]

#Data frame of stability parameters
data.nm <- data[191:375]
colnames(data.nm) <-1:185

#Create a subdataframe for only two mutants and wildtype PZAse
subdata.nm <- data.nm[c(1,9,21),]
subdata.nm <- as.data.frame(t(subdata.nm))
Position <- 1:185
subdata.nm <- cbind(subdata.nm,Position)

###########################
## Fluctuations profiles ##
###########################

#Create a new data set for plotting the fluctuation profile of three PZAses
subdata.nm <- gather(subdata.nm,key='Mutation',value='Fluctuation',-Position,factor_key = F)
g1 <- ggplot(subdata.nm,mapping = aes(x=Position,y=Fluctuation,color=Mutation)) + 
  geom_line(size=0.7) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.line = element_line(colour = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=seq(1,185,23), limits = c(1,185)) + 
  ylab(bquote('Fluctuation ('~ring(A)~')'))

################
## Loop 51-71 ##
################

# Creating a new variable that indicates if kcat is greater than 50
median.kcat <- round(median(data.k[,1]),2)
median.kcat <- 50
data.nm.k <- cbind(data.k[1],data.nm)
data.nm.g <- mutate(.data=data.nm.k,Kcat=as.factor(ifelse(data.nm.k[,1]>=median.kcat,
                                                          paste('\u2265 ',median.kcat,sep=''),
                                                          paste('< ',median.kcat,sep=''))))

#Function to plot fluctuation boxplots for each position 
plot_aa <- function(data,i=1,f=length(data)-1,d){
  data <- data[,c(1,(i+1):(f+1))]
  data <- gather(data,key='Position',value='Fluctuation',-Kcat,factor_key = T)
  ggplot(data,aes(x=Position,y=Fluctuation , fill = Kcat)) + 
    geom_boxplot(outlier.shape = 21,size=0.5) + theme(axis.title = element_text(size=12),
                                             axis.text = element_text(size=10),
                                             legend.text = element_text(size=12),
                                             legend.title = element_text(size=12))+
    scale_x_discrete(breaks=seq(i,f,d)) +
    scale_fill_discrete(name = bquote('Relative-'~k[cat])) +
    ylab(bquote('Fluctuation ('~ring(A)~')'))
}

#Plotting just the flap region
g2 <- plot_aa(data.nm.g,51,f=71,1)

#Save figure 1
grob <- arrangeGrob(g1, g2)
ggsave('fig1.tiff',grob,width = 6, height = 6, dpi = 300)

########################################################################
## Boxplot for the whole fluctuation profile for each group of PZAses ##
########################################################################

#Separate the dataset with kinetic information
data.nm.k.2 <- filter(data.nm.k, Kcat>=median.kcat)
data.nm.k.2 <- mutate(.data=data.nm.k.2,Kcat=as.factor(ifelse(data.nm.k.2[,1]>=median.kcat,
                                                              paste('\u2265 ',median.kcat,sep=''),
                                                              paste('< ',median.kcat,sep=''))))
data.nm.k.1 <- filter(data.nm.k, Kcat<median.kcat)
data.nm.k.1 <- mutate(.data=data.nm.k.1,Kcat=as.factor(ifelse(data.nm.k.1[,1]>=median.kcat,
                                                              paste('\u2265 ',median.kcat,sep=''),
                                                              paste('< ',median.kcat,sep=''))))
#Plot the mutants with kcat >= median
mnbp_s <- gather(data.nm.k.2,key='Position',value='Fluctuation',-Kcat,factor_key = T)
g3 <- ggplot(mnbp_s,aes(x=Position,y=Fluctuation, fill=Kcat)) + 
  geom_boxplot(outlier.shape = 21, fill='#00BFC4', size=0.3) + 
  theme(axis.title = element_text(size=12),axis.text = element_text(size=8)) +
  scale_x_discrete(breaks=seq(1,185,8)) + ylab(bquote('Fluctuation ('~ring(A)~')'))

#Plot the mutants with kcat < 25%
mnbp_r <- gather(data.nm.k.1,key='Position',value='Fluctuation',-Kcat,factor_key = T)
g4 <- ggplot(mnbp_r,aes(x=Position,y=Fluctuation, fill=Kcat)) + 
  geom_boxplot(outlier.shape = 21, fill='#F8766D', size=0.3) + 
  theme(axis.title = element_text(size=12), axis.text = element_text(size=8)) +
  scale_x_discrete(breaks=seq(1,185,8)) + ylab(bquote('Fluctuation ('~ring(A)~')'))

#Save supplementary figure 1
grob <- arrangeGrob(g3, g4)
ggsave('figsup1.tiff',grob,width = 6, height = 6, dpi = 300)

######################################################
## Perform a non parametric test in the flap region ##
######################################################

colnames(data.nm.g) <- c('Kcat',paste('x',1:185,sep=''))
position <- c()
p.values <- c()
conf.int.lower <- c()
conf.int.upper <- c()
for(i in colnames(data.nm.g[-1])){
  test <- wilcox.test(as.formula(paste(i,"~Kcat",sep="")),
                       alternative = 'two.sided',data.nm.g,conf.int=T)
  position <- c(position,i)
  conf.int.lower <- c(conf.int.lower,test$conf.int[1])
  conf.int.upper <- c(conf.int.upper,test$conf.int[2])
  p.values <- c(p.values,test$p.value)
  adj.p.values <- p.adjust(p.values, method ='BH')
}

test.table <-  cbind(position,conf.int.lower,conf.int.upper,p.values,adj.p.values)
write.table(test.table,'statistical_test_results.csv',sep=',',quote = F, row.names = F)

##################################
## Generalized statistical test ##
##################################

y <- c()
x <- c()
for(i in 1:(length(data.nm.k[,1])-5)){
  cutoff <- sort(data.nm.k[,1], decreasing = T)[i]
  data.nm.grouped <- mutate(.data=data.nm.k,Kcat=as.factor(ifelse(data.nm.k[,1]>=cutoff,'S','R')))
  colnames(data.nm.grouped) <- c('Kcat',paste('x',1:185,sep=''))
  p.values <- c()
  for(j in colnames(data.nm.grouped[-1])){
    p.values <- c(p.values,wilcox.test(as.formula(paste(j,"~Kcat",sep="")),
                        alternative = 'two.sided',data.nm.grouped)$p.value)
  }
  positions <- which(p.values < 0.05)
  y <- c(y,rep(table(data.nm.grouped[1])[1],length(positions)))
  x <- c(x,positions)
}

scale <- round(sort(data.nm.k[,1])[c(8,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)],2)
scale <- as.character(scale)
data.toplot = data.frame(N_resistant=y,position=x)
data.toplot[,1] <- as.factor(data.toplot[,1])
levels(data.toplot[,1]) <- 1:100
data.toplot[,1] <- as.integer(data.toplot[,1])
data.toplot[,2] <- as.factor(data.toplot[,2])
#levels(data.toplot[,2]) <- as.character(1:185)
ggplot(data.toplot,aes(x = N_resistant, y=position )) + geom_tile(color='black', fill='grey') +
  scale_x_continuous(breaks=1:26, labels = c(7,9,11:34),
                     sec.axis = dup_axis(labels = scale,name = '%kcat cutoff'))
#+  scale_y_continuous(1:185)

