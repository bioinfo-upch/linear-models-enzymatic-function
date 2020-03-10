#Molecular dynamics protein in water
echo -e "15"$ |gmx pdb2gmx -f *.pdb -o protein.gro -water spce
gmx editconf -f protein.gro -o protein_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top
gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr
echo -e "13"$ | gmx genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
#minimizacion
gmx grompp -f minim.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -nb gpu -deffnm em
#temperatura
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -nb gpu -deffnm nvt
#presion
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -nb gpu -deffnm npt
#dinamica
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -v -nb gpu -deffnm md_0_1

