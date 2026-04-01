import os
import re

# --- CONFIGURATION ---
half_pdb_folder = "/path/to/NAD(P)-binding/FirstHalf/folder/"
Ross_bab_folder = "/path/to/NAD(P)-binding/BAB/folder/"
output_folder = "/path/to/output/folder/"
log_file = os.path.join(output_folder, "topology_mismatch_log.txt")
faulty_list_path = "faulty_pdbs.txt" # The file containing the PDB codes to fix

os.makedirs(output_folder, exist_ok=True)

#The script is based on utilising the previously identified phosphate-binding BAB region and the NAD First Half to find the strand 2-3 segment
#The segment starts before the second beta strand of the phosphate-binding BAB and ends at the end of the NAD First Half


def get_pdb_code(filename): 
    #Extracts the 4-character PDB code from various naming styles
    match = re.search(r'([a-z0-9]{4})', filename.lower())
    return match.group(1) if match else None

def get_strand2_start(bab_path):
    #Finds the start of strand2 within the phosphate-binding BAB segments from NAD proteins
    strands = []
    with open(bab_path, 'r') as f:
        for line in f:
            if line.startswith("SHEET"):
                raw_val = line[21:27].strip()
                digits = re.sub(r'\D', '', raw_val)
                if digits:
                    strands.append(int(digits))
    strands.sort()
    return strands[1] if len(strands) >= 2 else None



# --- MAIN LOOP ---
half_files = [f for f in os.listdir(half_pdb_folder) if f.endswith(".pdb") and "_1_" in f] #Picks Half Files that are first in the structure "1"

for h_file in half_files:
    pdb_code = get_pdb_code(h_file)
    if not pdb_code:
        continue

    # Step 1: Find matching BAB file
    bab_file = next((f for f in os.listdir(Ross_bab_folder) if pdb_code in f.lower() and f.endswith(".pdb")), None)
    if not bab_file:
        print(f"Skipping {pdb_code}: No matching BAB file found.")
        continue

    # Step 2: Get start coordinate
    s2_start_res = get_strand2_start(os.path.join(Ross_bab_folder, bab_file))
    if s2_start_res is None: continue
    
    segment_start = s2_start_res - 1
    
    new_atoms = []
    new_ss = []
    
    # Step 3: Parse Half-PDB and extract segment
    with open(os.path.join(half_pdb_folder, h_file), 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_digits = re.sub(r'\D', '', line[21:27].strip())
                if res_digits and int(res_digits) >= segment_start:
                    new_atoms.append(line)
            elif line.startswith(("SHEET", "HELIX")):
                ss_digits = re.sub(r'\D', '', line[21:27].strip())
                if ss_digits and int(ss_digits) >= segment_start:
                    new_ss.append(line)

    # Write output
    output_name = f"{pdb_code}_BAB_strand2-3.pdb"
    with open(os.path.join(output_folder, output_name), 'w') as out_f:
        out_f.writelines(new_ss)
        out_f.writelines(new_atoms)
    
    
print(f"\nFinished!")
