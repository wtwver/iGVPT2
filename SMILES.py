import numpy as np
from rdkit.Chem import AddHs, MolFromSmiles, MolToXYZBlock, AllChem

def OrcaWrapper(theory, basis, *args, parallel=None):
    # header = f'! {" ".join([theory, basis, *args])} bohrs verytightscf ' + '{job_type:}\n'
    header = [
        '! RIJCOSX RI-B2PLYP D3BJ def2-tzvpp def2-tzvpp/C TIGHTSCF', 
         "! Opt NumFreq", 
        "* xyz 0 1" 
        ]
    mol = AddHs(MolFromSmiles("O"))
    AllChem.EmbedMolecule(mol)
    atom = MolToXYZBlock(mol).split("\n\n")[1]
    header.append(atom)
    header.append('*')

    return header

SMILESS = [
    "O"
#  "O=C(O)c1ccc(C(O)=O)cc1",
#  "C=C",
#  "ClC=C",
#  "C=CC",
# "c1ccccc1C=C",
# "CC(Cl)CC(Cl)CC(Cl)CC(Cl)CC(Cl)CC(Cl)"
]

orca = OrcaWrapper(theory='ri-mp2', basis='def2-tzvpp', parallel=4)

with open("test.inp","w") as f:
    for line in orca:
        print(line, file= f)

