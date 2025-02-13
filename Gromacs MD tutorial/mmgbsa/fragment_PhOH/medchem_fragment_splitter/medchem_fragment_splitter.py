import sys, getopt, os
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import traceback

def decomposition(argv):
    only_smarts = False
    picture = False
    write_fragments = False
    try:
        opts, args = getopt.getopt(argv, "l:w:s:pxa",
                                   ["ligand=", "workdir=", "smarts=", "picture", "write_frag", "auto"])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-l', '--ligand'):
            lig = arg
        elif opt in ('-w', '--workdir'):
            workdir = arg
        elif opt in ('-a', '--auto'):
            only_smarts =True
        elif opt in ('-p', '--picture'):
            picture =True
        elif opt in ('-x', '--write_frag'):
            write_fragments =True
        elif not only_smarts:
            if opt in ('-s', '--smarts'):
                smarts_path = arg
    smarts = list()
    if not only_smarts:
        with open(smarts_path, 'r') as smarts_reader:
            smarts_lines = smarts_reader.readlines()
        for smart_line in smarts_lines:
            smarts.append(smart_line.strip())

    try:
        if lig.endswith(".pdb"):
            molecule = Chem.MolFromPDBFile(lig)
        elif lig.endswith(".mol"):
            molecule = Chem.MolFromMol2File(lig)
        elif lig.endswith(".smi"):
            with open(lig, 'r') as smi_reader:
                smi = smi_reader.readline().strip().split(r'\s+')[0]
                molecule = Chem.MolFromSmiles(smi)
        elif lig.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(lig)
            for mol in suppl:
                molecule = mol
        else:
            molecule = Chem.MolFromSmiles(lig)

        atomPairs = []
        if not only_smarts:
            for index in range(len(smarts)):
                smartFragment = Chem.MolFromSmiles(smarts[index])
                atomPairsLoc = molecule.GetSubstructMatches(smartFragment)
                for atomPair in atomPairsLoc:
                    atomPairs.append(atomPair)
        else:
            for index in range(len(rules)):
                atomPairsLoc = molecule.GetSubstructMatches(Chem.MolFromSmarts(rules[index]))
                for nb in range(len(index_all[index])):
                    index1 = index_all[index][nb][0]
                    index2 = index_all[index][nb][1]
                    for atomPair in atomPairsLoc:
                        atomIndex1 = atomPair[index1]
                        atomIndex2 = atomPair[index2]
                        atomPairs.append((atomIndex1, atomIndex2))

        if (len(atomPairs) == 0):
            return
        atomPairs = list(set(atomPairs))
        bonds = list()
        if not only_smarts:
            flag = False
            for atomIndexes1 in atomPairs:
                for atomIndexes2 in atomPairs:
                    if atomIndexes1 == atomIndexes2:
                        continue
                    else:
                        for a1 in atomIndexes1:
                            for a2 in atomIndexes2:
                                if a1 == a2:
                                    flag = True
                                    if len(atomIndexes1) > len(atomIndexes2):
                                        if atomIndexes2 in atomPairs:
                                            atomPairs.remove(atomIndexes2)
                                    else:
                                        if atomIndexes1 in atomPairs:
                                            atomPairs.remove(atomIndexes1)
                                        break
                            if flag:
                                break
                    flag = False
        for atomIndexes in atomPairs:
            for a1 in atomIndexes:
                for a2 in atomIndexes:
                    if a1 == a2:
                        continue
                    else:
                         bond = molecule.GetBondBetweenAtoms(a1, a2)
                         if bond != None:
                             bonds.append(bond.GetIdx())
        bonds = list(set(bonds))
        if len(bonds)>0:
            cutMolecule = Chem.FragmentOnBonds(molecule, bonds, addDummies=False)
            cutMolecule.ClearComputedProps()
            for atom in cutMolecule.GetAtoms():
                if not atom.IsInRing():
                    if atom.GetIsAromatic():
                        atom.SetIsAromatic(False)
            try:
                tfragments = Chem.GetMolFrags(cutMolecule, asMols=True)
            except Chem.rdchem.AtomValenceException:
                cutMolecule = add_nitrogen_charges(cutMolecule)
                tfragments = Chem.GetMolFrags(cutMolecule, asMols=True)
            fragments = list()
            for f in tfragments:
                f = Chem.RemoveHs(f)
                f = Chem.AddHs(f, addCoords=True, addResidueInfo=True)
                fragments.append(f)
            fragment_count = 1
            if write_fragments:
                for fragment in fragments:
                    fragment_name = "ligand_" + str(fragment_count) + "_" + Descriptors.rdMolDescriptors.CalcMolFormula(
                        fragment)
                    Chem.MolToPDBFile(fragment, workdir + os.sep + fragment_name + ".pdb")
                    fragment_count += 1
            if picture:
                ifragments = Chem.GetMolFrags(cutMolecule, asMols=False)
                get_color_structure(molecule, fragments, ifragments, bonds, workdir+os.sep+os.path.basename(lig))
        else:
            get_color_structure(molecule, list(), list(), list(), workdir+os.sep+os.path.basename(lig))
    except:
        traceback.print_exc()
        print(lig)

def get_color_structure(molecule, fragments, ifragments, break_bonds, structure_output):
    AllChem.Compute2DCoords(molecule)
    imageStructure = rdMolDraw2D.MolDraw2DSVG(500, 250)
    imagePNGStructure = rdMolDraw2D.MolDraw2DCairo(500, 300)
    drawOpt = imageStructure.drawOptions()
    drawOpt.highlightBondWidthMultiplier = 16
    drawOpt.useBWAtomPalette()

    drawOptPNG = imagePNGStructure.drawOptions()
    drawOptPNG.highlightBondWidthMultiplier = 16
    drawOptPNG.useBWAtomPalette()
    bondsList = list()
    colorb = {}
    for bbond in break_bonds:
        bondsList.append(bbond)
        colorb[bbond] = (1.0, 0.0, 0.0)
    highlightAtoms = list()
    col_atoms = {}
    atom_radii = {}
    red_c = 0.6
    green_c = 0.84
    blue_c = 0.6
    current_color = (red_c, green_c, blue_c)

    for index in range(len(fragments)):
        highlightAtoms.extend(ifragments[index])
        numAtomsInFragment = len(ifragments[index])
        for atom1 in range(numAtomsInFragment):
            col_atoms[ifragments[index][atom1]] = current_color
            atom_radii[ifragments[index][atom1]] = 0.45
            for atom2 in range(atom1 + 1, numAtomsInFragment):
                try:
                    bond = molecule.GetBondBetweenAtoms(ifragments[index][atom1], ifragments[index][atom2])
                    if (bond != None):
                        bond_id = bond.GetIdx()
                        colorb[bond_id] = (red_c, green_c, blue_c)
                        bondsList.append(bond_id)
                except Exception:
                    pass
    
    # draw all at once
    rdMolDraw2D.PrepareAndDrawMolecule(imageStructure, molecule, highlightAtoms=highlightAtoms,
                                           highlightAtomRadii=atom_radii,
                                           highlightBonds=bondsList,
                                           highlightAtomColors=col_atoms, highlightBondColors=colorb)
    rdMolDraw2D.PrepareAndDrawMolecule(imagePNGStructure, molecule, highlightAtoms=highlightAtoms,
                                           highlightAtomRadii=atom_radii,
                                           highlightBonds=bondsList,
                                           highlightAtomColors=col_atoms, highlightBondColors=colorb)
    imageStructure.FinishDrawing()
    imagePNGStructure.FinishDrawing()
    svg = imageStructure.GetDrawingText().replace('svg:', '')
    print(svg, file=open(structure_output+".svg", 'w'))
    imagePNGStructure.WriteDrawingText(structure_output+".png")

def add_nitrogen_charges(m):
    m.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(m)
    if not ps:
        Chem.SanitizeMol(m)
        return m
    for p in ps:
        if p.GetType()=='AtomValenceException':
            at = m.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                at.SetFormalCharge(1)
    Chem.SanitizeMol(m)
    return m

rules = [
    "[R]-!@[$([CX4;H2,H1,H0])]",  # 1
    "[a]-!@[$([NX3;H1,H0]),$([OX2;H0]),$([SX2;H0])]-!@[$([C;H2,H1,H0]);!$([CX3]=[OX1])]",  # 2
    "[a]-!@[$([NX3;H1,H0]),$([OX2;H0]),$([SX2;H0])]-!@[a]",  # 3
    "[a]-!@[$([CX3]=[OX1,NX2,SX1,CX3])]-!@[$([CX4;H2,H1,H0])]",  # 4
    "[c]-!@[$([CX3]=[OX1,NX2,SX1,$([CX3;H2])])]-!@[c]",  # 5.1
    "[n]-!@[$([CX3]=[OX1,NX2,SX1,$([CX3;H2])])]-!@[c]",  # 5.2
    "[$([CX4;H2,H1,H0])]-!@[CX3](=[OX1])[OX2;H0]",  # 6
    "[$([CX4;H2,H1,H0])]-!@[OX2;H0][CX3](=[OX1])",  # 7
    "[a]-!@[CX3](=[OX1])O-!@[$([CX4;H2,H1,H0])]",  # 8
    "[a]-!@[CX3](=[OX1])O-!@[a]",  # 9
    "[a]-!@[NX2;H0]=[NX2;H0]-!@[$([CX4;H2,H1,H0])]",  # 10
    "[a]-!@[NX2;H0]=[NX2;H0]-!@[a]",  # 11
    "[a]-!@[NX3;H1]-!@[$([CX3;H0](=[OX1]))]-!@[$([CX4;H2,H1,H0])]",  # 12
    "[a]-!@[$([CX3;H0](=[OX1]))]-!@[NX3;H1]-!@[$([CX4;H2,H1,H0])]",  # 13
    "[a]-!@[NX3;H1]-!@[$([CX3;H0](=[OX1]))]-!@[a]",  # 14
    "[a]-!@[$([CX3;H0](=[OX1]))]-!@[NX3;H1]-!@[a]",  # 15
    "[a]-!@[SX4](=[OX1])(=[OX1])[NX3;H1]-!@[$([CX4;H2,H1,H0])]",  # 16
    "[a]-!@[NX3;H1][SX4](=[OX1])(=[OX1])-!@[$([CX4;H2,H1,H0])]",  # 17
    "[a]-!@[SX4](=[OX1])(=[OX1])[NX3;H1]-!@[a]",  # 18
    "[a]-!@[NX3;H1][SX4](=[OX1])(=[OX1])-!@[NX3;H1]-!@[a]",  # 19
    "[$([CX4;H2,H1,H0])]-!@[NX3][CX3](=[OX1])",  # 20
    "[$([CX4;H2,H1,H0])]-!@[CX3](=[OX1])[NX3]",  # 21
    "[$([CX4;H2,H1,H0])]-!@[$([NX3;H1,H0])]",  # 22
    "[$([CX4;H2,H1,H0])]-!@[$([OX2;H0])]",  # 23
    "[$([CX4;H2,H1,H0])]-!@[$([SX2;H0])]",  # 24
    "[$([CX4;H2,H1,H0])]-!@[SX4](=[OX1])(=[OX1])[NX3;H1]",  # 25
    "[$([CX4;H2,H1,H0])]-!@[NX3;H1][SX4](=[OX1])(=[OX1])"  # 26
]

index_all = [
    [(0, 1)],  # 1
    [(1, 2)],  # 2
    [(0, 1), (1, 2)],  # 3
    [(1, 2)],  # 4
    [(0, 1), (1, 2)],  # 5.1
    [(1, 2)],  # 5.2
    [(0, 1)],  # 6
    [(0, 1)],  # 7
    [(3, 4)],  # 8
    [(0, 1), (3, 4)],  # 9
    [(2, 3)],  # 10
    [(0, 1), (2, 3)],  # 11
    [(2, 3)],  # 12
    [(2, 3)],  # 13
    [(0, 1), (2, 3)],  # 14
    [(0, 1), (2, 3)],  # 15
    [(4, 5)],  # 16
    [(4, 5)],  # 17
    [(0, 1), (4, 5)],  # 18
    [(0, 1), (4, 5)],  # 19
    [(0, 1)],  # 20
    [(0, 1)],  # 21
    [(0, 1)],  # 22
    [(0, 1)],  # 23
    [(0, 1)],  # 24
    [(0, 1)],  # 25
    [(0, 1)],  # 26
]

if __name__ == '__main__':
    decomposition(sys.argv[1:])
