import os
from Bio.PDB import PDBParser
import math

#Function to calculate the angle between the vector from PA to N9A and PA to N10
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


#Function to process each pdb file in a folder and extract positions of FAD if present
def process_fad_folder(folder_path, output_file):
    parser = PDBParser(QUIET=True)
    results = []

    for fname in os.listdir(folder_path):
        if not fname.endswith(".pdb"):
            continue
        pdb_path = os.path.join(folder_path, fname)
        structure = parser.get_structure(fname, pdb_path)
        found = False

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname().strip() == "FAD":
                        atoms = {"N9A": None, "PA": None, "N10": None}   #Extract position of 3 atoms N9A, PA and N10
                        for atom in residue:
                            name = atom.get_name().strip().replace('"', '')
                            if name in atoms:
                                atoms[name] = atom.get_coord()
                        if all(value is not None for value in atoms.values()):
                            angle = calculate_angle(atoms["N9A"], atoms["PA"], atoms["N10"])
                            results.append((fname, angle))
                            found = True
                            break
                if found:
                    break
            if found:
                break
        if not found:
            results.append((fname, "N/A"))

    with open(output_file, "w") as f:
        for fname, angle in results:
            f.write(f"{fname}\t{angle}\n")

# === Run for your FAD folder ===
process_fad_folder("/path/to/FAD/Folder/", "/path/to/output/file.txt")
