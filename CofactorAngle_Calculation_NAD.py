import os
from Bio.PDB import PDBParser
import math

#Function to calculate the angle between the vector from PA to N9A and PA to C4N
def calculate_angle(a, b, c):
    def vector(p1, p2):
        return [p2[i] - p1[i] for i in range(3)]
    def dot(v1, v2):
        return sum(v1[i]*v2[i] for i in range(3))
    def norm(v):
        return math.sqrt(sum(coord**2 for coord in v))
    
    ba = vector(b, a)
    bc = vector(b, c)
    cosine_angle = dot(ba, bc) / (norm(ba) * norm(bc))
    angle_rad = math.acos(cosine_angle)
    return math.degrees(angle_rad)


#Function to process each pdb file in a folder and extract positions of NAD if present
def process_nad_folder(folder_path, output_file):
    parser = PDBParser(QUIET=True)
    results = []

    for fname in os.listdir(folder_path):
        if not fname.endswith(".pdb"):
            continue
        pdb_path = os.path.join(folder_path, fname)
        structure = parser.get_structure(fname, pdb_path)
        found = False

        atoms = {"N9A": None, "PA": None, "C4N": None}  #Extract position of 3 atoms N9A, PA and C4N
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname().strip() == "NAD":
                        for atom in residue:
                            name = atom.get_name().strip().replace('"', '')
                            if name in atoms:
                                atoms[name] = atom.get_coord()

        if all(v is not None for v in atoms.values()):
            angle = calculate_angle(atoms["N9A"], atoms["PA"], atoms["C4N"])
            results.append((fname, angle))
        else:
            results.append((fname, "N/A"))


    with open(output_file, "w") as f:
        for fname, angle in results:
            f.write(f"{fname}\t{angle}\n")

# === Run for your NAD folder ===
process_nad_folder("/path/to/NAD/Folder", "/path/to/ouput/file.txt")
