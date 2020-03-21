 
# Authors: Eduardo Gushiken & Rydberg Supo Escalante

####################
## Load libraries ##
####################

packages <- c('dplyr', 'purrr')
for(i in packages){library(package=i,character.only = T)}

###############################
## Read PDBs and coordinates ##
###############################

# Name of the PDBs
mutants = c("WT","A102V","A171T","A46V","C14G","D12A","D12G","D136G","D49N","F58L","F94L",
            "G24D","G78C","H51R","H57R","H71Y","K48T","L116P","L172P","L4S","M175V","P54L",
            "P62L","P62R","Q10P","R29P","T135P","T142A","T160K","T76P","V125F","V139A",
            "V180F","W119L","Y34D","Y64D")

# Name of positions
variables = as.character(1:185)

# List to store the coordinates for each structure (alpha carbon)
coordinates <- list(coordinates_x = data.frame(),coordinates_y = data.frame(),
                    coordinates_z = data.frame())

# List to store the resultant vector of each c
resultant_vector <- list(resultant_vector_x = data.frame(), resultant_vector_y = data.frame(),
                         resultant_vector_z = data.frame())

# This loop iterates on each PDB and extracts coordinates and resultant vector
for (structure in mutants){
  
  # Store alpha carbons 
  c_alpha <- list(c1=as.numeric(),c2=as.numeric(),c3=as.numeric())
  
  # Store atoms different to alpha carbons
  c_no_alpha <- list(c1_no_a = list(),c2_no_a = list(),c3_no_a = list())
  
  # Open pdb
  structure = paste(structure,".pdb",sep="")
  con = file(structure, open="r")
  
  # Number of amino-acid
  aa_old = 0
  
  # Read pdb
  for (i in readLines(con)){
    
    # Extract information from pdb
    x = substr(i,1,4) #ATOM
    y = substr(i,14,15) #CA 
    aa = as.numeric(substr(i,24,26)) #Amino-acid number
    cxyz=as.list(as.numeric(c(substr(i,31,38),substr(i,39,46),substr(i,47,54)))) #Coordinates
    
    # If alpha carbon, store coordinates in c_alpha
    if (y == "CA" & x == "ATOM"){ 
      c_alpha <- map2(c_alpha,cxyz,c)
      aa_old = aa
    }
    
    # If not alpha carbon, store coordinates in c_no_alpha
    if (y != "CA" & x == "ATOM"){ 
      if (aa != aa_old){
        c_no_alpha <- map2(c_no_alpha, cxyz, append)
      } else {
        c_no_alpha <- pmap(list(c_no_alpha, list(aa),list(c), cxyz), modify_at)
      }
      aa_old = aa
    }
    
  }
  
  # Close pdb
  close(con)
  
  # Calculate resultant vector for each position
  resultant_vector <- map2(resultant_vector, 
                           modify(pmap(list(map(c_no_alpha,map_dbl,sum), 
                                            map(c_no_alpha,map_dbl,length), c_alpha), 
                                       function(x,y,z){x-y*z}),~as.data.frame(t(.x))),rbind)

  # Store alpha carbon coordinates
  coordinates <- map2(coordinates,modify(c_alpha,~as.data.frame(t(.x))),rbind)
}

#Output: Resultant vector & Ca coordinates

###########################
## Calculate Baricenters ##
###########################

# Baricenter AS (D8, K96, C138)
baricenter_3_as = list(baricenter_x = data.frame(), baricenter_y = data.frame(), 
                    baricenter_z = data.frame())
ba_sa <- map(map(coordinates,select,8,96,138),~rowSums(.x)/3)
baricenter_3_as <- map2(baricenter_3_as,map(ba_sa,cbind),rbind)

# Baricenter MCS (D49, H51, H71)
baricenter_3_mcs = list(baricenter_x = data.frame(), baricenter_y = data.frame(), 
                       baricenter_z = data.frame())
ba_mcs <- map(map(coordinates,select,49,51,71),~rowSums(.x)/3)
baricenter_3_mcs <- map2(baricenter_3_mcs,map(ba_mcs,cbind),rbind)

# Baricenter MCS (D49, H51, H57, H71) - T point
baricenter_4_mcs = list(baricenter_t_x = data.frame(), baricenter_t_y = data.frame(), 
                       baricenter_t_z = data.frame())
ba_mcs_t <- map(map(coordinates,select,49,51,57,71),~rowSums(.x)/4)
baricenter_4_mcs <- map2(baricenter_4_mcs,map(ba_mcs_t,cbind),rbind)

# Output: baricenter_3_as, baricenter_3_mcs, baricenter_4_mcs

###############################################
## Calculate the Normal Vector of AS and MCS ##
###############################################

# List to store the coordinates of the normal vector
Normal_Vector <- list(normal_x = data.frame(), normal_y = data.frame(), normal_z = data.frame())

# Calculate normal vector of the AS
vm_as <- map2(map(coordinates, select, 96, 138), map(map(coordinates, select, 8),as_vector), ~.x-.y)
vm_as <- map(map2(list(vm_as$coordinates_y,vm_as$coordinates_z,vm_as$coordinates_x),
         list(vm_as$coordinates_z,vm_as$coordinates_x,vm_as$coordinates_y),
         function(a,b){a[,2]*b[,1]-b[,2]*a[,1]}),cbind)

# Calculate normal vector of the MCS (D49, H51, H71)
vm_mcs <- map2(map(coordinates, select, 51, 71), map(map(coordinates, select, 49),as_vector), ~.x-.y)
vm_mcs <-  map(map2(list(vm_mcs$coordinates_y,vm_mcs$coordinates_z,vm_mcs$coordinates_x),
                    list(vm_mcs$coordinates_z,vm_mcs$coordinates_x,vm_mcs$coordinates_y),
                    function(a,b){a[,2]*b[,1]-b[,2]*a[,1]}),cbind)

# Store normal vectors
Normal_Vector <- map2(Normal_Vector,map2(vm_as,vm_mcs,cbind),rbind)

# Calculate the constant D of the plane equation 
d <- map(coordinates,select,8,49)
d <- map2(Normal_Vector,d,~.x*.y)
Normal_Vector$normal_d <- pmap_dfc(list(list(d[[1]]),list(d[[2]]),list(d[[3]])),
                                   function(x,y,z){-(x+y+z)})

# Name the columns of each data frame
Normal_Vector <- map(Normal_Vector,function(x){colnames(x) = c("D_AS","D_MCS");x})

# Output: Normal_Vector

###############################################
## Calculate I and P points for each residue ##
###############################################

# List to store the coordinates of the point I of the AS
I_as <- list(I_as_x = data.frame(), I_as_y = data.frame(), I_as_z = data.frame())

# List to store the coordinates of the point I of the MCS
I_mcs <- list(I_mcs_x = data.frame(), I_mcs_y = data.frame(), I_mcs_z = data.frame())

# List to store the coordinates of the point P of the AS
P_as <- list(P_as_x = data.frame(), P_as_y = data.frame(), P_as_z = data.frame())

# List to store the coordinates of the point P of the MCS
P_mcs <- list(P_mcs_x = data.frame(), P_mcs_y = data.frame(), P_mcs_z = data.frame())

# Normal vector and resultant vector of each position will be used to calculate this points

# Normal vector of the AS
as_n <- map(map(Normal_Vector,select,1),as_vector)

# Normal vector of the MCS
mcs_n <- map(map(Normal_Vector,select,2),as_vector)

# Calculation of I_as
as_t <- map2(as_n[1:3],coordinates,function(x,y){x*y})
as_t <- pmap_dfc(list(list(as_t[[1]]), list(as_t[[2]]), list(as_t[[3]])), function(x,y,z){x+y+z})
as_t_part1 <- -(as_t + as_n[[4]])
as_t <- map2(as_n[1:3],resultant_vector,function(x,y){x*y})
as_t_part2 <- pmap_dfc(list(list(as_t[[1]]), list(as_t[[2]]), list(as_t[[3]])),function(x,y,z){x+y+z})
as_t <- as_t_part1/as_t_part2
I_as <- map2(coordinates,map(resultant_vector, function(n){as_t*n}),~.x+.y)

# Calculation of I_mcs
mcs_t <- map2(mcs_n[1:3],coordinates,function(x,y){x*y})
mcs_t <- pmap_dfc(list(list(mcs_t[[1]]), list(mcs_t[[2]]), list(mcs_t[[3]])), function(x,y,z){x+y+z})
mcs_t_part1 <- -(mcs_t + mcs_n[[4]])
mcs_t <- map2(mcs_n[1:3],resultant_vector,function(x,y){x*y})
mcs_t_part2 <- pmap_dfc(list(list(mcs_t[[1]]), list(mcs_t[[2]]), list(mcs_t[[3]])),function(x,y,z){x+y+z})
mcs_t <- mcs_t_part1/mcs_t_part2
I_mcs <- map2(coordinates,map(resultant_vector, function(n){mcs_t*n}),~.x+.y)

# Calculation of P_as
as_s <- map2(as_n[1:3],as_n[1:3],~.x*.y)
as_s <- pmap(list(list(as_s[[1]]),list(as_s[[2]]),list(as_s[[3]])), function(x,y,z){x+y+z})
as_s <- as_t_part1/as_s
P_as <- map2(coordinates,map(as_n[1:3],function(n){n*as_s}),~.x+.y)

# Calculation of P_mcs
mcs_s <- map2(mcs_n[1:3],mcs_n[1:3],~.x*.y)
mcs_s <- pmap(list(list(mcs_s[[1]]),list(mcs_s[[2]]),list(mcs_s[[3]])), function(x,y,z){x+y+z})
mcs_s <- mcs_t_part1/mcs_s
P_mcs <- map2(coordinates,map(mcs_n[1:3],function(n){n*mcs_s}),~.x+.y)

# Outputs:  I_as, I_mcs, P_as, P_mcs.

################################
## Function to write Datasets ##
################################

# Inputs:
# list1 and list2: Lists containing three dataframes, one for each coordinate
# mutants: Rownames of the output dataframe
# variables: Position names
# name: tag to include in the name of the output

distance <- function(list1, list2, mutants, variables, name){
  
  # Prepare datasets
  list1 <- ifelse(sapply(list1, length)==185, list1, map(list1,as_vector))
  list2 <- ifelse(sapply(list2, length)==185, list2, map(list2,as_vector))
  
  # Calculate distance
  diff <- map2(list1,list2, function(x1,x2){(x2-x1)^2})
  dist <- as.data.frame(pmap_dfc(list(list(diff[[1]]),list(diff[[2]]),list(diff[[3]])),
                                 function(x,y,z){sqrt(x+y+z)}))
  # Name rows and columns
  rownames(dist) <- mutants
  
  if(length(dist) == 185){
    colnames(dist) <- variables
  } else {
    colnames(dist) <- ('D')
  }
  
  # Write dataset
  name = paste("distance_",name,".csv",sep="")
  write.csv(dist,file=name)
  return(dist)
}

# Outputs: A csv file with the distances for each structure and a dataframe with the same information

####################
## Write Datasets ##
####################

# Use the distace function to calculate distances involving the lists baricenter_3_as, 
# baricenter_3_mcs, baricenter_4_mcs, coordinates, I_mcs, I_as, P_mcs, P_as
D8 <- distance(coordinates,map(coordinates,select,8),mutants,variables,'to_aa_8')
D49 <- distance(coordinates,map(coordinates,select,49),mutants,variables,'to_aa_49')
H51 <- distance(coordinates,map(coordinates,select,51),mutants,variables,'to_aa_51')
H57 <- distance(coordinates,map(coordinates,select,57),mutants,variables,'to_aa_57')
H71 <- distance(coordinates,map(coordinates,select,71),mutants,variables,'to_aa_71')
K96 <- distance(coordinates,map(coordinates,select,96),mutants,variables,'to_aa_96')
C138 <- distance(coordinates,map(coordinates,select,138),mutants,variables,'to_aa_138')
B_C_as <- distance(coordinates,baricenter_3_as,mutants,variables,'B_C_as')
B_C_mcs <- distance(coordinates,baricenter_3_mcs,mutants,variables,'B_C_mcs')
B_I_as  <- distance(I_as, baricenter_3_as,mutants,variables,"B_I_as")
B_I_mcs <- distance(I_mcs,baricenter_3_mcs,mutants,variables,"B_I_mcs")
B_P_as  <- distance(baricenter_3_as, P_as,mutants,variables,"B_P_as")
B_P_mcs <- distance(baricenter_3_mcs, P_mcs,mutants,variables,"B_P_mcs")
I_C_as  <- distance(coordinates,I_as,mutants,variables,"I_C_as")
I_C_mcs <- distance(coordinates,I_mcs,mutants,variables,"I_C_mcs")
I_P_as  <- distance(P_as,I_as, mutants,variables,"I_P_as")
I_P_mcs <- distance(P_mcs,I_mcs,mutants,variables,"I_P_mcs")
P_C_as  <- distance(P_as,coordinates,mutants,variables,"P_C_as")
P_C_mcs <- distance(P_mcs,coordinates,mutants,variables,"P_C_mcs")
C_T_mcs <- distance(coordinates,baricenter_4_mcs,mutants,variables,"C_T_mcs")
I_T_mcs <- distance(I_mcs,baricenter_4_mcs, mutants,variables,"I_T_mcs")
P_T_mcs <- distance(P_mcs,baricenter_4_mcs, mutants,variables,"P_T_mcs")
Ias_Imcs <- distance(I_mcs, I_as, mutants,variables,"Ias_Imcs")
Pas_Pmcs <- distance(P_mcs,P_as,mutants,variables,"Pas_Pmcs")

# Create a combined dataset gathering all the information
data <- cbind(D8,D49,H51,H57,H71,K96,C138,B_C_as,B_C_mcs,B_I_as,B_I_mcs,B_P_as,B_P_mcs,I_C_as,
              I_C_mcs,I_P_as,I_P_mcs,P_C_as,P_C_mcs,C_T_mcs,I_T_mcs,P_T_mcs,Ias_Imcs,Pas_Pmcs)
colnames(data) <- paste('x',1:length(data),sep='')
write.csv(data,file='geometrical_descriptors.csv')

