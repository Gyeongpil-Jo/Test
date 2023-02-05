#!/bin/bash

# requirements:
# amber tools & program_1.nab file
source /apps/jejoong/amber20/amber.sh
# gromacs
module load gromacs/2020.2
# /apps/jejoong/jejoong_utils/*.pl files

# The only input
seq=cgccgccgcgggttcctctaggagattctcac
nseq=32
echo "Insert ...$seq..."


# 1. PDB file of DNA using amber tools.
sed "s/REPLACE/$seq/" ../program_1.nab > program_1.nab

nab -o program_1 program_1.nab 
./program_1


# 2. cut DNA that we want (remove first and last nucleotides)
gmx make_ndx -f dna.pdb -o tmp.ndx <<+
del 1-30
ri 2-$(($nseq+1))
ri $(($nseq+4))-$((2*$nseq+4-1))
1 | 2
q
+

# a PDB file with DNA structure of the given sequence with 10-bp periodicity.

# write chain 1
gmx editconf -f dna.pdb -o dna1.pdb -n tmp.ndx -resnr 1 <<+
1
+

# write chain 2
gmx editconf -f dna.pdb -o dna2.pdb -n tmp.ndx -resnr 1 <<+
2
+

#DT B   1

# combine  two chains into one PDB file (dna.pdb)
perl -i -p -e 's/D([ATGC]) (.)   1/D$1J$2   1/' dna?.pdb
perl -i -p -e 's/D([ATGC])3/D$1 /' dna?.pdb
cat dna[12].pdb | /bin/grep -E "^ATOM" > dna.pdb
echo TER >> dna.pdb


# 3. Prepare Gromacs files using dna.pdb
gmx pdb2gmx -f dna.pdb   -merge all -ignh  -chainsep id -o dna0.pdb -p dna <<+
14
1
+

# cut top 25 lines and bottom 25 lines in dna.itp
N0=$(wc -l dna.top | awk '{print $1}')
N1=$(  echo "$N0-25" | bc -l)
N2=$(  echo "$N1-27" | bc -l)
head -n $N1 dna.top | tail -n $N2 > dna.itp

amber_dna_pbc_hbonds.pl  dna.itp # output will be dna.pbc.itp

# moleculetype is complete

echo 0 | gmx editconf -f dna0.pdb -o tmp1.pdb -princ -box 6  6  $(echo "$nseq*0.338" | bc -l)
gmx editconf -f tmp1.pdb  -o tmp2.pdb -c -center 0 0 0
gmx editconf -f tmp2.pdb  -o tmp3.pdb -rotate 0 90 0

# now, dna is in parallel to z and centered at origin.

gmx editconf -f tmp3.pdb  -o tmp4.pdb -c 

# twist by 360 degree over contour
#phi=$(echo "-360/($nseq*3.4)" | bc -l)
pdb2twist4pbc.pl --phi=-0.7 tmp4.pdb > dna.pdb

cadnano2pdb2enm.pl dna.pdb > dna.enm.itp

perl -i -p -e 's/90.00  90.00  90.00/90.00  90.00  60.00/' dna.pdb
gmx editconf -f dna.pdb -o dna.pdb -c

exit















#####################
# place arg5 along minor groove
#for i in $(seq.pl 0 9)
#do
#    gmx editconf -f ../lys.pdb -o lys$i.pdb -translate 0 0 $(echo "3*$i*0.34" | bc -l)
#    gmx editconf -f lys$i.pdb -o lys$i.pdb -rotate 0 0 $(echo "(360*3/32)*3*$i" | bc -l)
#    gmx editconf -f lys$i.pdb -o lys$i.pdb -translate 4.6 2.5 0
#done
#awk '/^ATOM/' lys?.pdb  > lys.pdb
awk '/^CRYS/' dna.pdb > tmp6.pdb
awk '/^ATOM/' dna.pdb ../arg5.pdb >> tmp6.pdb
gmx solvate -cp tmp6.pdb -o tmp7.pdb -cs

#rm -f lys?.pdb *#*

NNARG=$(grep "CA  ARG" tmp7.pdb | wc -l)
NNSOL=$(grep OW tmp7.pdb | wc -l)
NNNA=$(grep NA tmp7.pdb | wc -l)
NNCL=$(grep CL tmp7.pdb | wc -l)

rm -f topol.top
cp ../top.sample topol.top
perl -i -p -e "s/NNARG/$NNARG/" topol.top
perl -i -p -e "s/NNSOL/$NNSOL/" topol.top
perl -i -p -e "s/NNNA/$NNNA/" topol.top
perl -i -p -e "s/NNCL/$NNCL/" topol.top

gmx grompp -f ../mini.mdp -c tmp7.pdb -maxwarn 40

gmx genion -np  79 -nn 20 -o conf.pdb <<+
SOL
+

NNSOL=$(grep OW conf.pdb | wc -l)
NNNA=$(grep NA conf.pdb | wc -l)
NNCL=$(grep CL conf.pdb | wc -l)

rm -f topol.top
cp ../top.sample topol.top
perl -i -p -e "s/NNARG/$NNARG/" topol.top
perl -i -p -e "s/NNSOL/$NNSOL/" topol.top
perl -i -p -e "s/NNNA/$NNNA/" topol.top
perl -i -p -e "s/NNCL/$NNCL/" topol.top

gmx grompp -f ../mini.mdp -maxwarn 45

gmx mdrun -nt 4 -c conf.pdb

../prep.sh

exit





gmx make_ndx -f topol.tpr <<+
a P
a NZ & ri 66
a NZ & ri 69
a NZ & ri 72
a NZ & ri 75
a NZ & ri 78
name 23 NZ1
name 24 NZ2
name 25 NZ3
name 26 NZ4
name 27 NZ5
r D* & !a H*
name 33 DNA-H
q
+

gmx grompp -f ../grompp.0.2fs.mdp -n -o topol.0.2fs.tpr -maxwarn 45
#grompp -f ../grompp.mdp -n -o topol.tpr -maxwarn 5
gmx grompp -f ../pull.mdp -n -o topol.tpr -maxwarn 45

../now.sh

rm -f *#* *tmp* *step* lys?.pdb

