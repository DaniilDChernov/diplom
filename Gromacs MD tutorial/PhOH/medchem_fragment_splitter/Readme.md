MedChem fragment splitter (MedChemFrag)
=========================
Split a drug like molecule into fragments meaningful from the viewpoint of Medicinal Chemistry.

**Description**
The following options are available in script:
- -l - path to ligand file (available molecule format: pdb, mol2, smi, sdf);
- -r - path to receptor file, option if needed to save a complex with each separate fragment in a one file;
- -w - path to working directory;
- -o - prefix for output file name;
- -f - output molecule format for fragments (available: pdb, mol, smi, sdf)
- -s - path to file with SMARTS expression;
- -p - no parameter, option if needed a picture;
- -x - no parameter, option if needed to save fragment-files
- -a - no parameter, option if needed to be decomposed by inner rules
- -d - no parameter, option if needed to save file with information about which atom in the original molecule is in which fragment. The file is a mol file with additional "V" flags ("V {atom index} {number of fragment}")

The application can work in several modes. The first mode decomposes molecules by the built-in set of SMARTS-expressions. It is achieved using parameter ‘-s’ and a file with defined SMARTS-expressions:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -s ..SMARTS-expressions -o prefix_output_name -f output_format -x -p

What is more, these expressions must describe non-overlaid fragments. The second mode lets to decompose molecules by the rules described in current article. It is achieved using parameter ‘-a’:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -o prefix_output_name -f output_format -x -p -a

The parameter ‘-r’ is used to generate a complex of a molecule's with receptor:

>./defragmentation.py -l ..path-to-ligand -r ..path-to-receptor -w ..path-to-work-directory -o prefix_output_name -f output_format -p -a

The parameter ‘-x’ is used to write files of fragments:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -o prefix_output_name -f output_format -x -a

The parameter ‘-p’ is used to generate a picture of a molecule's structure with colored fragments:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -o prefix_output_name -f output_format -p -a

The parameter ‘-d’ is used to generate a mol-file of molecule with information about which atom in the original molecule is in which fragment:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -o prefix_output_name -f output_format -x -a -d

The output file is a mol-file with additional "V" flags ("V {atom index} {number of fragment}"). In order to extract this information using rdkit, you need to use the "GetProp" method:

> num_fragment = mol.GetAtomWithIdx(atom_index).GetProp("molFileValue")

To produce the output at least one parameter ‘-x’ or ‘-p’ should be used.

The SDF format can be used to decompose multiple molecules. If needed to write fragments to one file, then use sdf as the output file format.

**Installation and requirements**
- Python3
- RDKit version 2020.03 or later

An Ubuntu 22.04 based Docker container is supplied and is run for the test by issuing in the base directory of the project

>./bin/docker_test.sh

**Examples of application**
The test supplied splits two known drugs with DrugBank IDs DB00176 and DB08865 into fragments.
The pictures and the fragments could be found in files directory after the run.

**Scientific references**
https://doi.org/10.3390/biophysica3020024

