MedChem fragment splitter (MedChemFrag)
=========================
Split a drug like molecule into fragments meaningful from the viewpoint of Medicinal Chemistry.

**Descriptopn**
The following options are available in script:
-l - path to ligand file;
-w - path to working directory;
-s - path to file with SMARTS expression;
-p - no parameter, option if needed a picture;
-x - no parameter, option if needed to save fragment-files
-a - no parameter, option if needed to be decomposed by inner rules

The application can work in several modes. The first mode decomposes molecules by the built-in set of SMARTS-expressions. It is achieved using parameter ‘-s’ and a file with defined SMARTS-expressions:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -s ..SMARTS-expressions -x -p

What is more, these expressions must describe non-overlaid fragments. The second mode lets to decompose molecules by the rules described in current article. It is achieved using parameter ‘-a’:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -x -p -a

The parameter ‘-x’ is used to write files of fragments:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -x -a

The parameter ‘-p’ is used to generate a picture of a molecule's structure with colored fragments:

>./defragmentation.py -l ..path-to-ligand -w ..path-to-work-directory -p -a

To produce the output at least one parameter ‘-x’ or ‘-p’ should be used.

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

