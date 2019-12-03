#!/bin/bash

# Author: Basilio Cieza

############################################
# Eigenvalue format change
############################################

head -n 46  eigenvalues_modif.txt > temp.value.1.txt
tail -n 40  temp.value.1.txt > temp.value.2.txt
awk '{print $2}' temp.value.2.txt > temp.value.3.txt
mv temp.value.3.txt eigenvalues_modif.txt

############################################
# Eigenvector format change
############################################

tail -n -1 eigenvectors_modif.TXT > temp.vector.1.txt
mv  temp.vector.1.txt  eigenvectors_modif.TXT
rm temp*
mv eigenvectors_modif.TXT eigenvector.txt
mv eigenvalues_modif.txt eigenvalue.txt
rm centers_modif.txt  HESSIAN_modif.txt beta_modif.txt 

############################################

Rscript FORMAT.R

############################################

############################################
# Changing format to be read by fluctuation script
############################################

#change format for cordinate
sed 's/"//g' vector.output.txt > temp0.txt
awk '{print $2}' temp0.txt > temp1.txt
awk '{print $3}' temp0.txt > temp2.txt
awk '{print $4}' temp0.txt > temp3.txt
paste temp1.txt temp2.txt temp3.txt > temp4.txt
sed 1d temp4.txt >> temp5.txt

#Insert column with atom data
sed 's/"/" /g' vector.output.txt > col.temp1.txt
sed 1d col.temp1.txt >> col.temp2.txt
awk '{print $1}' col.temp2.txt > col.temp3.txt
sed 's/"/   2 ASP   7 CA /g' col.temp3.txt > col.temp4.txt

#Jjoin atom data with coordinate data
paste col.temp4.txt temp5.txt > output.ready.txt
rm *temp*
rm flucxyz_modif.txt vector.output.txt eigenvector.txt
python editing.py 
rm output.ready.txt
