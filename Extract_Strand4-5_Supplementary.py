import os
import re

# --- CONFIGURATION ---
sec_half_folder = "/path/to/NAD(P)-binding/SecondHalf/Rossmann/PDBs/" 
output_folder = "/path/to/output/folder/"

os.makedirs(output_folder, exist_ok=True)

#The extraction is based on utilising the previously identified NAD(P)-binding Second half of all PDBs to find the start of the Strand4-5 segment

def get_pdb_code(filename):
    #Extracts the 4-character PDB code from various naming styles
    match = re.search(r'([a-z0-9]{4})', filename.lower())
    return match.group(1) if match else None

def get_strand5_end(pdb_path):
    #Finds the end residue of the 2nd strand record (Strand 5) found in the Second Half file.
    strand_ends = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("SHEET"):
                # PDB columns for end residue: 34-37
                try:
                    res_str = line[33:37].strip()
                    # Clean insertion codes if present
                    res_num = int(''.join(filter(str.isdigit, res_str)))
                    strand_ends.append(res_num)
                except ValueError:
                    continue
    
    # In the second half file, the 1st strand found is Strand 4, 
    # the 2nd strand found is Strand 5.
    return strand_ends[1] if len(strand_ends) >= 2 else None

# --- MAIN LOOP ---
sec_files = [f for f in os.listdir(sec_half_folder) if f.endswith(".pdb") and "_1_" in f] #Picks Half Files that are first in the structure "1"


for s_file in sec_files:
    pdb_code = get_pdb_code(s_file)
    if not pdb_code: continue
    
    path = os.path.join(sec_half_folder, s_file)
    s5_end_res = get_strand5_end(path)
    
    if s5_end_res is None:
        print(f"Skipping {pdb_code}: Could not find a 5th strand in the second half.")
        continue

    # Define the breakpoint: 1 residue after the end of strand 5
    segment_end = s5_end_res + 1
    
    new_atoms = []
    new_ss = []
    
    with open(path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    res_num = int(''.join(filter(str.isdigit, line[22:26].strip())))
                    if res_num <= segment_end:
                        new_atoms.append(line)
                except ValueError:
                    continue
            elif line.startswith(("SHEET", "HELIX")):
                try:
                    # Include SS record if its start is before our endpoint
                    ss_start = int(''.join(filter(str.isdigit, line[22:27].strip())))
                    if ss_start <= segment_end:
                        new_ss.append(line)
                except ValueError:
                    continue

    # Write output
    output_name = f"{pdb_code}_BAB_strand4-5.pdb"
    with open(os.path.join(output_folder, output_name), 'w') as out_f:
        out_f.writelines(new_ss)
        out_f.writelines(new_atoms)

print(f"Extraction of Strands 4-5 complete. Files saved to {output_folder}")


