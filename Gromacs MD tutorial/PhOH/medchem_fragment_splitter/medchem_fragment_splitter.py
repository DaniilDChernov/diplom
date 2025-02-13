import sys, getopt, os, json
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import traceback

def decomposition(argv):
    only_smarts = False
    picture = False
    write_fragments = False
    receptor = None
    rec = None
    produce_fragment_mol = False

    # ligands - path to ligand file
    # workdir - path to working directory (for save fragments files)
    # smarts - path to file with SMARTS
    # output_prefix - prefix of output file name
    # output_format - output format
    # receptor - path to receptor file
    # pictures - if to need save fragments image
    # write_frag - if to need save fragments files
    # auto - if to need use internal smarts dictionary
    try:
        opts, args = getopt.getopt(argv, "l:w:s:o:f:r:pxad",
              ["ligand=", "workdir=", "smarts=", "output_prefix=", "output_format=", "receptor=", 
              "picture", "write_frag", "auto", "produce_fragment_mol"])
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
        elif opt in ('-o', '--output_prefix'):
            output_prefix = arg
        elif opt in ('-f', '--output_format'):
            output_format = arg
        elif opt in ('-r', '--receptor'):
            receptor = arg
        elif opt in ('-s', '--smarts'):
            smarts_path = arg
        elif opt in ('-d', '--produce_fragment_mol'):
            produce_fragment_mol = True
    smarts = list()
    if not only_smarts:
        with open(smarts_path, 'r') as smarts_reader:
            smarts_lines = smarts_reader.readlines()
        for smart_line in smarts_lines:
            smarts.append(smart_line.strip())
    else:
        with open(os.path.join(os.path.dirname(__file__), 'rules.json'), 'r') as json_reader:
            rules_json = json.load(json_reader)
            rules = rules_json["smarts"]
            index_all = rules_json["indexes"]

    try:
        if receptor!=None:
            rec = Chem.MolFromPDBFile(receptor)
        if not lig.endswith(".sdf"):
            if lig.endswith(".pdb"):
                molecule = Chem.MolFromPDBFile(lig)
                #molecule = Chem.MolFromPDBFile(lig, removeHs=False)
            elif lig.endswith(".mol"):
                #molecule = Chem.MolFromMolFile(lig)
                molecule = Chem.MolFromMolFile(lig, removeHs=False)
            elif lig.endswith(".mol2"):
                #molecule = Chem.MolFromMol2File(lig)
                molecule = Chem.MolFromMol2File(lig, removeHs=False)
            elif lig.endswith(".smi"):
                with open(lig, 'r') as smi_reader:
                    smi = smi_reader.readline().strip().split(r'\s+')[0]
                    molecule = Chem.MolFromSmiles(smi)
            else:
                molecule = Chem.MolFromSmiles(lig)
            if receptor!=None and not molecule.GetConformer().Is3D():
                print("Ligand has not 3D coordinates. The complex will not be recorded!")
                rec=None
            #processing_molecules(lig, only_smarts, smarts, rules, index_all, molecule, write_fragments, picture, output_format, output_prefix, workdir, rec)
            # shulga: produce fragments using explicit Hs and produce pictures using implicit Hs - they are much prettier
            processing_molecules(lig, only_smarts, smarts, rules, index_all, molecule, write_fragments, 
                False, # picture 
                output_format, output_prefix, workdir, rec, produce_fragment_mol)
            if picture:
                pic_molecule = Chem.RemoveHs(molecule)
                processing_molecules(lig, only_smarts, smarts, rules, index_all, 
                pic_molecule, #molecule, 
                False, #write_fragments=False, 
                picture, output_format, output_prefix, workdir, rec,
                False #produce_fragment_mol
                )
        else:
            suppl = Chem.SDMolSupplier(lig, removeHs=False)
            fragments_sdf = list()
            count = 1
            for molecule in suppl:
                #fragments_molecule = processing_molecules(lig, only_smarts, smarts, rules, index_all, molecule, write_fragments, picture,
                #                     output_format, output_prefix+"_"+str(count), workdir, None)
                fragments_molecule = processing_molecules(lig, only_smarts, smarts, rules, index_all, molecule, write_fragments, 
                    False, #picture,
                    output_format, output_prefix+"_"+str(count), workdir, None, produce_fragment_mol)
                if picture:
                    pic_molecule = Chem.RemoveHs(molecule)
                    processing_molecules(lig, only_smarts, smarts, rules, index_all, pic_molecule, 
                        False, # write_fragments, 
                        picture,
                        output_format, output_prefix+"_"+str(count), workdir, None, 
                        False #produce_fragment_mol
                        )

                count += 1
                if output_format=="sdf":
                    fragments_sdf.extend(fragments_molecule)
        if lig.endswith(".sdf") and output_format=="sdf":
            with Chem.SDWriter(output_prefix+"."+output_format) as sd_writer:
                for f in fragments_sdf:
                    sd_writer.write(f)

    except:
        traceback.print_exc()
        print(lig)



def processing_molecules(lig, only_smarts, smarts, rules, index_all, molecule, write_fragments, picture, output_format, output_prefix, workdir, rec, produce_fragment_mol):
    try:
        molecule2=Chem.Mol(molecule)
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
        if len(bonds) > 0:
            cutMolecule = Chem.FragmentOnBonds(molecule, bonds, addDummies=False)
            #cutMolecule.ClearComputedProps()
            #Chem.SanitizeMol(cutMolecule)
            #for atom in cutMolecule.GetAtoms():
            #    if not atom.IsInRing():
            #        if atom.GetIsAromatic():
            #            atom.SetIsAromatic(False)
            try:
                #tfragments = Chem.GetMolFrags(cutMolecule, asMols=True)
                tfragments = Chem.GetMolFrags(cutMolecule, asMols=True, sanitizeFrags=True)
            except Chem.rdchem.AtomValenceException:
                #print("here exception")
                cutMolecule = add_nitrogen_charges(cutMolecule)
                tfragments = Chem.GetMolFrags(cutMolecule, asMols=True)
            fragments = list()
            for f in tfragments:
                # replace radicals with implicit hydrogens
                for a in f.GetAtoms():
                    re = a.GetNumRadicalElectrons()
                    if re != 0:
                        a.SetNumRadicalElectrons(0)
                        a.SetNumExplicitHs(a.GetNumExplicitHs() + re)
                        a.UpdatePropertyCache()
                #f = Chem.RemoveHs(f)
                #f = Chem.RemoveHs(f, implicitOnly=True)
                #Chem.SanitizeMol(f, catchErrors=True)
                f = Chem.AddHs(f, addCoords=True, addResidueInfo=True)
                #f = Chem.AddHs(f, addCoords=True, addResidueInfo=True, explicitOnly=True)
                fragments.append(f)
#                # DEBUG
#                print(f"fragment {Chem.MolToSmiles(f)}")
#                for  atom  in  f.GetAtoms():
#                    print(atom.GetAtomicNum(),  atom.GetNumImplicitHs(),
#                        atom.GetNumExplicitHs(),  atom.GetNumRadicalElectrons(), atom.GetFormalCharge())
            if picture:
                ifragments = Chem.GetMolFrags(cutMolecule, asMols=False)
                get_color_structure(molecule, ifragments, bonds, workdir + os.sep + output_prefix)

            if produce_fragment_mol:
                # DEBUG
                # print("here")
                ifragments = Chem.GetMolFrags(cutMolecule, asMols=False)
                numfragcount=1
                for ifg in ifragments:
                    for i in ifg:
                        molecule2.GetAtomWithIdx(i).SetProp("molFileValue", str(numfragcount))
                    numfragcount+=1
                Chem.MolToMolFile(molecule2, workdir + os.sep + output_prefix+"_out.mol")

            if write_fragments:
                fragment_count = 1
                sdf_mols = None
                fragments_part_sdf = list()
                if output_format == "sdf" and not lig.endswith(".sdf"):
                    sdf_mols = Chem.SDWriter(workdir + os.sep + output_prefix + ".sdf")
                for fragment in fragments:
                    fragment_name = output_prefix + "_" + str(
                        fragment_count)
                    if rec != None:
                        fragment = Chem.CombineMols(rec, fragment)
                    if output_format == "pdb":
                        Chem.MolToPDBFile(fragment, workdir + os.sep + fragment_name + ".pdb")
                    elif output_format == "smi":
                        with open(fragment_name, 'w') as fragment_writer:
                            fragment_writer.write(Chem.MolToSmiles(fragment))
                    elif output_format == "mol":
                        if rec!=None:
                            fragment = Chem.CombineMols(rec, fragment)
                        Chem.MolToMolFile(fragment, workdir + os.sep + fragment_name + ".mol")
                    elif output_format == "sdf":
                        if lig.endswith(".sdf"):
                            fragment.SetProp('_Name', fragment_name)
                            fragments_part_sdf.append(fragment)
                        else:
                            sdf_mols.write(fragment)
                    fragment_count += 1
                if lig.endswith(".sdf") and output_format=="sdf":
                    return fragments_part_sdf
        else:
            get_color_structure(molecule, list(), list(), workdir + os.sep + output_prefix)
    except:
        traceback.print_exc()
        print(lig)

def get_color_structure(molecule, ifragments, break_bonds, structure_output):
    ## shulga
    #molecule = Chem.RemoveHs(molecule)
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

    red_c = 0.6
    green_c = 0.84
    blue_c = 0.6

    col_atoms = {}
    current_color = (red_c, green_c, blue_c)
    atom_radii = {}

    for index in range(len(ifragments)):
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
    atom_indexes = list()
    for frag_index in ifragments:
        atom_indexes.extend(frag_index)
    rdMolDraw2D.PrepareAndDrawMolecule(imageStructure, molecule, highlightAtoms=atom_indexes,
                                       highlightAtomRadii=atom_radii,
                                       highlightBonds=bondsList,
                                       highlightAtomColors=col_atoms, highlightBondColors=colorb)
    rdMolDraw2D.PrepareAndDrawMolecule(imagePNGStructure, molecule, highlightAtoms=atom_indexes,
                                       highlightAtomRadii=atom_radii,
                                       highlightBonds=bondsList,
                                       highlightAtomColors=col_atoms, highlightBondColors=colorb)
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

if __name__ == '__main__':
    decomposition(sys.argv[1:])
