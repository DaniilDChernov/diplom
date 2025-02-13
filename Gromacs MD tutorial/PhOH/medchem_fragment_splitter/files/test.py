import os, difflib, sys

no_clean = False
if len(sys.argv) > 1:
    if sys.argv[1] == "--no_clean":
        no_clean = True

### TEST DB08865 ###

os.system("python3 ../medchem_fragment_splitter.py -l DB08865.smi -w . -f pdb -o DB08865 -p -x -a -d")
pic_DB08865_png = os.path.exists("DB08865.png")
pic_DB08865_svg = os.path.exists("DB08865.svg")
if pic_DB08865_png:
    print("Draw PNG Test 1 OK!")
else:
    print("Draw PNG Test 1 ERROR!")
if pic_DB08865_svg:
    print("Draw SVG Test 2 OK!")
else:
    print("Draw SVG Test 2 ERROR!")

fragment_1_DB08865 = "DB08865_1.pdb"
fragment_2_DB08865 = "DB08865_2.pdb"
fragment_3_DB08865 = "DB08865_3.pdb"
fragment_4_DB08865 = "DB08865_4.pdb"
mol_DB08854="DB08865_out.mol"

fragment_1_DB08865_out = "DB08865_1_out.pdb"
fragment_2_DB08865_out = "DB08865_2_out.pdb"
fragment_3_DB08865_out = "DB08865_3_out.pdb"
fragment_4_DB08865_out = "DB08865_4_out.pdb"
mol_DB08854_out="DB08865_out_out.mol"

diff_fragment_1_DB08865 = difflib.SequenceMatcher(None, open(fragment_1_DB08865).readlines(), open(fragment_1_DB08865_out).readlines()).ratio()
diff_fragment_2_DB08865 = difflib.SequenceMatcher(None, open(fragment_2_DB08865).readlines(), open(fragment_2_DB08865_out).readlines()).ratio()
diff_fragment_3_DB08865 = difflib.SequenceMatcher(None, open(fragment_3_DB08865).readlines(), open(fragment_3_DB08865_out).readlines()).ratio()
diff_fragment_4_DB08865 = difflib.SequenceMatcher(None, open(fragment_4_DB08865).readlines(), open(fragment_4_DB08865_out).readlines()).ratio()

diff_mol_DB08865 = difflib.SequenceMatcher(None, open(mol_DB08854).readlines(), open(mol_DB08854_out).readlines()).ratio()

if (diff_fragment_1_DB08865+diff_fragment_2_DB08865+diff_fragment_3_DB08865+diff_fragment_4_DB08865+diff_mol_DB08865)==5.0:
    print("Decomposition Test 3 OK!")
else:
    print("Decomposition Test 3 ERROR!")

if not no_clean:
    os.remove("DB08865.png")
    os.remove("DB08865.svg")
    os.remove(fragment_1_DB08865)
    os.remove(fragment_2_DB08865)
    os.remove(fragment_3_DB08865)
    os.remove(fragment_4_DB08865)
    os.remove(mol_DB08854)


### TEST DB00176 ###
os.system("python3 ../medchem_fragment_splitter.py -l DB00176.pdb -w . -f mol -o DB00176 -p -x -a -d")
pic_DB00176_png = os.path.exists("DB00176.png")
pic_DB00176_svg = os.path.exists("DB00176.svg")
if pic_DB00176_png:
    print("Draw PNG Test 4 OK!")
else:
    print("Draw PNG Test 4 ERROR!")
if pic_DB00176_svg:
    print("Draw SVG Test 5 OK!")
else:
    print("Draw SVG Test 5 ERROR!")

fragment_1_DB00176 = "DB00176_1.mol"
fragment_2_DB00176 = "DB00176_2.mol"
fragment_3_DB00176 = "DB00176_3.mol"
fragment_4_DB00176 = "DB00176_4.mol"
fragment_5_DB00176 = "DB00176_5.mol"
mol_DB00176="DB00176_out.mol"

fragment_1_DB00176_out = "DB00176_1_out.mol"
fragment_2_DB00176_out = "DB00176_2_out.mol"
fragment_3_DB00176_out = "DB00176_3_out.mol"
fragment_4_DB00176_out = "DB00176_4_out.mol"
fragment_5_DB00176_out = "DB00176_5_out.mol"
mol_DB00176_out="DB00176_out_out.mol"

diff_fragment_1_DB00176 = difflib.SequenceMatcher(None, open(fragment_1_DB00176).readlines(), open(fragment_1_DB00176_out).readlines()).ratio()
diff_fragment_2_DB00176 = difflib.SequenceMatcher(None, open(fragment_2_DB00176).readlines(), open(fragment_2_DB00176_out).readlines()).ratio()
diff_fragment_3_DB00176 = difflib.SequenceMatcher(None, open(fragment_3_DB00176).readlines(), open(fragment_3_DB00176_out).readlines()).ratio()
diff_fragment_4_DB00176 = difflib.SequenceMatcher(None, open(fragment_4_DB00176).readlines(), open(fragment_4_DB00176_out).readlines()).ratio()
diff_fragment_5_DB00176 = difflib.SequenceMatcher(None, open(fragment_5_DB00176).readlines(), open(fragment_5_DB00176_out).readlines()).ratio()

diff_mol_DB00176 = difflib.SequenceMatcher(None, open(mol_DB00176).readlines(), open(mol_DB00176_out).readlines()).ratio()

if (diff_fragment_1_DB00176+diff_fragment_2_DB00176+diff_fragment_3_DB00176+diff_fragment_4_DB00176+diff_fragment_5_DB00176+diff_mol_DB00176)==6.0:
    print("Decomposition Test 6 OK!")
else:
    print("Decomposition Test 6 ERROR!")

if not no_clean:
    os.remove("DB00176.png")
    os.remove("DB00176.svg")
    os.remove(fragment_1_DB00176)
    os.remove(fragment_2_DB00176)
    os.remove(fragment_3_DB00176)
    os.remove(fragment_4_DB00176)
    os.remove(fragment_5_DB00176)
    os.remove(mol_DB00176)

os.system("python3 ../medchem_fragment_splitter.py -l DB00176.mol2 -w . -f sdf -o DB00176 -p -x -a -d")
pic_DB00176_png = os.path.exists("DB00176.png")
pic_DB00176_svg = os.path.exists("DB00176.svg")
if pic_DB00176_png:
    print("Draw PNG Test 7 OK!")
else:
    print("Draw PNG Test 7 ERROR!")
if pic_DB00176_svg:
    print("Draw SVG Test 8 OK!")
else:
    print("Draw SVG Test 8 ERROR!")

fragments_DB00176 = "DB00176.sdf"
mol_DB00176="DB00176_out.mol"

fragments_DB00176_out = "DB00176_out.sdf"
mol_DB00176_out="DB00176_mol2_out.mol"

diff_fragments_DB00176 = difflib.SequenceMatcher(None, open(fragments_DB00176).readlines(), open(fragments_DB00176_out).readlines()).ratio()
diff_mol_DB00176 = difflib.SequenceMatcher(None, open(mol_DB00176).readlines(), open(mol_DB00176_out).readlines()).ratio()

if (diff_fragments_DB00176+diff_mol_DB00176)==2.0:
    print("Decomposition Test 9 OK!")
else:
    print("Decomposition Test 9 ERROR!")

if not no_clean:
    os.remove("DB00176.png")
    os.remove("DB00176.svg")
    os.remove(fragments_DB00176)
    os.remove(mol_DB00176)

### TEST 3HF7 ###
#os.system("python3 ../medchem_fragment_splitter.py -l 3hf7_ligand.pdb -r 3hf7_receptor.pdb -w . -f pdb -o 3hf7 -p -x -a")
os.system("python3 ../medchem_fragment_splitter.py -l 3hf7_ligand.mol -r 3hf7_receptor.pdb -w . -f pdb -o 3hf7 -p -x -a -d")
pic_3hf7_png = os.path.exists("3hf7.png")
pic_3hf7_svg = os.path.exists("3hf7.svg")
if pic_3hf7_png:
    print("Draw PNG Test 10 OK!")
else:
    print("Draw PNG Test 10 ERROR!")
if pic_3hf7_svg:
    print("Draw SVG Test 11 OK!")
else:
    print("Draw SVG Test 11 ERROR!")

fragment_1_3hf7 = "3hf7_1.pdb"
fragment_2_3hf7 = "3hf7_2.pdb"
fragment_3_3hf7 = "3hf7_3.pdb"
fragment_4_3hf7 = "3hf7_4.pdb"
mol_3hf7="3hf7_out.mol"

fragment_1_3hf7_out = "3hf7_1_out.pdb"
fragment_2_3hf7_out = "3hf7_2_out.pdb"
fragment_3_3hf7_out = "3hf7_3_out.pdb"
fragment_4_3hf7_out = "3hf7_4_out.pdb"
mol_3hf7_out="3hf7_out_out.mol"

diff_fragment_1_3hf7 = difflib.SequenceMatcher(None, open(fragment_1_3hf7).readlines(), open(fragment_1_3hf7_out).readlines()).ratio()
diff_fragment_2_3hf7 = difflib.SequenceMatcher(None, open(fragment_2_3hf7).readlines(), open(fragment_2_3hf7_out).readlines()).ratio()
diff_fragment_3_3hf7 = difflib.SequenceMatcher(None, open(fragment_3_3hf7).readlines(), open(fragment_3_3hf7_out).readlines()).ratio()
diff_fragment_4_3hf7 = difflib.SequenceMatcher(None, open(fragment_4_3hf7).readlines(), open(fragment_4_3hf7_out).readlines()).ratio()

diff_mol_3hf7 = difflib.SequenceMatcher(None, open(mol_3hf7).readlines(), open(mol_3hf7_out).readlines()).ratio()

if (diff_fragment_1_3hf7+diff_fragment_2_3hf7+diff_fragment_3_3hf7+diff_fragment_4_3hf7+diff_mol_3hf7)==5.0:
    print("Decomposition Test 12 OK!")
else:
    print("Decomposition Test 12 ERROR!")

if not no_clean:
    os.remove("3hf7.png")
    os.remove("3hf7.svg")
    os.remove(fragment_1_3hf7)
    os.remove(fragment_2_3hf7)
    os.remove(fragment_3_3hf7)
    os.remove(fragment_4_3hf7)
    os.remove(mol_3hf7)

### Test DB ###
os.system("python3 ../medchem_fragment_splitter.py -l DB.sdf -w . -f sdf -o DBF -p -x -a -d")
pic_DBF_1_png = os.path.exists("DBF_1.png")
pic_DBF_1_svg = os.path.exists("DBF_1.svg")
pic_DBF_2_png = os.path.exists("DBF_2.png")
pic_DBF_2_svg = os.path.exists("DBF_2.svg")
if pic_DBF_1_png & pic_DBF_2_png:
    print("Draw PNG Test 13 OK!")
else:
    print("Draw PNG Test 13 ERROR!")
if pic_DBF_1_svg & pic_DBF_2_svg:
    print("Draw SVG Test 14 OK!")
else:
    print("Draw SVG Test 14 ERROR!")

fragments_DBF = "DBF.sdf"
mol_DBF_1 = "DBF_1_out.mol"
mol_DBF_2 = "DBF_2_out.mol"

fragments_DBF_out = "DBF_out.sdf"
mol_DBF_1_out = "DBF_1_out_out.mol"
mol_DBF_2_out = "DBF_2_out_out.mol"

diff_fragments_DBF = difflib.SequenceMatcher(None, open(fragments_DBF).readlines(), open(fragments_DBF_out).readlines()).ratio()
diff_mol_DBF_1 = difflib.SequenceMatcher(None, open(mol_DBF_1).readlines(), open(mol_DBF_1_out).readlines()).ratio()
diff_mol_DBF_2 = difflib.SequenceMatcher(None, open(mol_DBF_2).readlines(), open(mol_DBF_2_out).readlines()).ratio()

if (diff_fragments_DBF+diff_mol_DBF_1+diff_mol_DBF_2)==3.0:
    print("Decomposition Test 15 OK!")
else:
    print("Decomposition Test 15 ERROR!")

if not no_clean:
    os.remove("DBF_1.png")
    os.remove("DBF_2.png")
    os.remove("DBF_1.svg")
    os.remove("DBF_2.svg")
    os.remove("DBF.sdf")
    os.remove("DBF_1_out.mol")
    os.remove("DBF_2_out.mol")

