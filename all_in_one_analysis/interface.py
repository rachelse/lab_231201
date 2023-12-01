from Bio.PDB import *

def get_structure(name, path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, path)
    return structure

def get_interface(structure, threshold=8):
    c = list(structure[0].child_dict.keys())
    chain1 = structure[0][c[0]]
    chain2 = structure[0][c[1]]
    
    interface_pairs = []
    interface1 = [] 
    interface2 = []

    for r1 in chain1:
        for r2 in chain2:
            try:
                if r1['CA'] - r2['CA'] <= threshold:
                    interface_pairs.append((r1.get_id()[1], r2.get_id()[1]))
                    interface1.append(r1.get_id()[1])
                    interface2.append(r2.get_id()[1])
            except:
                continue

    interface1 = sorted(set(interface1))
    interface2 = sorted(set(interface2))

    interface = {'pairs' : interface_pairs, c[0] : interface1, c[1] : interface2}
    return interface

if __name__ == "__main__":
    path = "/home/seamustard52/dpi/benchmark/pdb/036679_039621.pdb"
    structure = get_structure("068352_068435", path)
    interface = get_interface(structure)
    print(interface)