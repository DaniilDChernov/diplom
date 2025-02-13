#!/bin/bash

python3 ../medchem_fragment_splitter.py -l DB08865.smi -w . -p -x -a
if [ -f DB08865.smi.svg ]; then
	if [ -s DB08865.smi.svg ]; then
		echo "Draw SVG Test 1 OK!"
	else
		echo "Draw SVG Test 1 Error!"
	fi
else
	echo "Draw SVG Test 1 Error!"
fi
if [ -f DB08865.smi.png ]; then
	if [ -s DB08865.smi.png ]; then
		echo "Draw PNG Test 2 OK!"
	else
		echo "Draw PNG Test 2 Error!"
	fi
else
	echo "Draw PNG Test 2 Error!"
fi

PDB_1_LIGAND=$(diff ligand_1_C2H4.pdb output_ligand_1_C2H4.pdb | wc -c)
PDB_2_LIGAND=$(diff ligand_2_C8H8N4O.pdb output_ligand_2_C8H8N4O.pdb | wc -c)
PDB_3_LIGAND=$(diff ligand_3_C5H11N.pdb output_ligand_3_C5H11N.pdb | wc -c)
PDB_4_LIGAND=$(diff ligand_4_C6H3Cl2F.pdb output_ligand_4_C6H3Cl2F.pdb | wc -c)
PDB_SUM=0
((PDB_SUM+=PDB_1_LIGAND))
((PDB_SUM+=PDB_2_LIGAND))
((PDB_SUM+=PDB_3_LIGAND))
((PDB_SUM+=PDB_4_LIGAND))

if [[ ${PDB_SUM} -eq '0' ]]; then
	echo "Decomposition Test 3 OK!"
else
	echo "Decomposition Test 3 Error!"
fi

# clean up
rm DB08865.smi.svg
rm DB08865.smi.png
rm ligand_1_C2H4.pdb
rm ligand_2_C8H8N4O.pdb
rm ligand_3_C5H11N.pdb
rm ligand_4_C6H3Cl2F.pdb

# Test 4-6
python3 ../medchem_fragment_splitter.py -l DB00176.smi -w . -p -x -a
if [ -f DB00176.smi.svg ]; then
	if [ -s DB00176.smi.svg ]; then
		echo "Draw SVG Test 1 OK!"
	else
		echo "Draw SVG Test 1 Error!"
	fi
else
	echo "Draw SVG Test 1 Error!"
fi

if [ -f DB00176.smi.png ]; then
	if [ -s DB00176.smi.png ]; then
		echo "Draw PNG Test 2 OK!"
	else
		echo "Draw PNG Test 2 Error!"
	fi
else
	echo "Draw PNG Test 2 Error!"
fi

PDB_1_LIGAND=$(diff ligand_1_CH4O.pdb output_ligand_1_CH4O.pdb | wc -c)
PDB_2_LIGAND=$(diff ligand_2_C4H10.pdb output_ligand_2_C4H10.pdb | wc -c)
PDB_3_LIGAND=$(diff ligand_3_C7H7NO.pdb output_ligand_3_C7H7NO.pdb | wc -c)
PDB_4_LIGAND=$(diff ligand_4_C2H7N.pdb output_ligand_4_C2H7N.pdb | wc -c)
PDB_5_LIGAND=$(diff ligand_5_CHF3.pdb output_ligand_5_CHF3.pdb | wc -c)
PDB_SUM=0
((PDB_SUM+=PDB_1_LIGAND))
((PDB_SUM+=PDB_2_LIGAND))
((PDB_SUM+=PDB_3_LIGAND))
((PDB_SUM+=PDB_4_LIGAND))
((PDB_SUM+=PDB_5_LIGAND))

if [[ ${PDB_SUM} -eq '0' ]]; then
	echo "Decomposition Test 6 OK!"
else
	echo "Decomposition Test 6 Error!"
fi

# clean up
rm DB00176.smi.svg
rm DB00176.smi.png
rm ligand_1_CH4O.pdb
rm ligand_2_C4H10.pdb
rm ligand_3_C7H7NO.pdb
rm ligand_4_C2H7N.pdb
rm ligand_5_CHF3.pdb
