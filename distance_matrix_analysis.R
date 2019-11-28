
# Author: Rydberg Supo Escalante

########################
## Charging libraries ##
########################

packages <- c('dplyr', 'caret', 'olsrr', 'car', 'rgl','reshape2','ggplot2','biobroom','tidyr',
              'dendextend','matrixStats','cluster','gridExtra','grid','lattice')
for(i in packages){library(package=i,character.only = T)}

#################################################
## Charging the datasets and distance analysis ##
#################################################

#Define jet colormap
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#Loading data
data_8 <- read.table('distance_to_aa_8.csv',sep=',',header = T, row.names = 1)[1,]
rownames(data_8) <- 'D8'
data_96 <- read.table('distance_to_aa_96.csv', sep = ',', header = T, row.names = 1)[1,]
rownames(data_96) <- 'K96'
data_138 <- read.table('distance_to_aa_138.csv', sep=',', header = T, row.names = 1)[1,]
rownames(data_138) <- 'C138'
data_49 <- read.table('distance_to_aa_49.csv', sep=',', header = T, row.names = 1)[1,]
rownames(data_49) <- 'D49'
data_51 <- read.table('distance_to_aa_51.csv', sep=',', header = T, row.names = 1)[1,]
rownames(data_51) <- 'H51'
data_57 <- read.table('distance_to_aa_57.csv', sep=',', header = T, row.names = 1)[1,]
rownames(data_57) <- 'H57'
data_71 <- read.table('distance_to_aa_71.csv', sep=',', header = T, row.names = 1)[1,]
rownames(data_71) <- 'H71'

#Create working data
data <- rbind(data_8,data_96,data_138,data_49,data_51,data_57,data_71)
data <- as.data.frame(t(data))
data <- round(data,2)
residues <- readLines('residues.txt')
rownames(data) <- residues
selected_residues <- readLines('selected_residues.txt')

#Generating the plots
D8 <- data[order(data[1],decreasing = T),]
D8 <- as.matrix(D8[1])
D8 <- melt(D8)
D8 <- D8[-c(1:165),]
h <- ifelse(D8[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
D8 <- ggplot(D8, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(legend.position='none',axis.title.x = element_blank(), 
        axis.title.y = element_text(size=30),axis.text = element_text(size=20),
        axis.text.y = element_text(color=h)) + geom_text(aes(label=value), size = 6) +
  ylab('Positions')

K96 <- data[order(data[2],decreasing = T),]
K96 <- as.matrix(K96[2])
K96 <- melt(K96)
K96 <- K96[-c(1:165),]
h <- ifelse(K96[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
K96 <- ggplot(K96, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(legend.position='none',axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text = element_text(size=20),  axis.text.y = element_text(color=h)) +
  geom_text(aes(label=value), size = 6) 

C138 <- data[order(data[3],decreasing = T),]
C138 <- as.matrix(C138[3])
C138 <- melt(C138)
C138 <- C138[-c(1:165),]
h <- ifelse(C138[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
C138 <- ggplot(C138, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(legend.position='none',axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text = element_text(size=20), axis.text.y = element_text(color=h))+
  geom_text(aes(label=value), size = 6) 

D49 <- data[order(data[4],decreasing = T),]
D49 <- as.matrix(D49[4])
D49 <- melt(D49)
D49 <- D49[-c(1:165),]
h <- ifelse(D49[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
D49 <- ggplot(D49, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(legend.position='none',axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text = element_text(size=20), axis.text.y = element_text(color=h)) +
  geom_text(aes(label=value), size = 6)

H51 <- data[order(data[5],decreasing = T),]
H51 <- as.matrix(H51[5])
H51 <- melt(H51)
H51 <- H51[-c(1:165),]
h <- ifelse(H51[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
H51 <- ggplot(H51, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(legend.position='none',axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text = element_text(size=20), axis.text.y = element_text(color=h)) +
  geom_text(aes(label=value), size = 6)

H57 <- data[order(data[6],decreasing = T),]
H57 <- as.matrix(H57[6])
H57 <- melt(H57)
H57 <- H57[-c(1:165),]
h <- ifelse(H57[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
H57 <- ggplot(H57, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(legend.position='none',axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text = element_text(size=20), axis.text.y = element_text(color=h)) +
  geom_text(aes(label=value), size = 6)

H71 <- data[order(data[7],decreasing = T),]
H71 <- as.matrix(H71[7])
H71 <- melt(H71)
H71 <- H71[-c(1:165),]
h <- ifelse(H71[,'Var1'] %in% selected_residues,yes = 'blue',no = 'black')
H71 <- ggplot(H71, aes(x=Var2, y = Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(7), limits = c(0,11)) +
  theme(axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text = element_text(size=20), axis.text.y = element_text(color=h),
        legend.text = element_text(size=20),legend.title = element_text(size=25),
        legend.key.size = unit(1, "cm")) + geom_text(aes(label=value), size = 6) +
  labs(fill = "Distance")

#Save the results
lay <- rbind(c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7))
grob <- arrangeGrob(D8, K96, C138, D49, H51, H57, H71, layout_matrix = lay,
                    top = textGrob('Distance matrix',gp=gpar(fontsize=30,font=3)))
ggsave('distance_heatmap.tiff',grob,width = 16, height = 16, dpi = 300)

