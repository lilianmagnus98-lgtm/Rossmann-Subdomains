import os
import re
from collections import Counter


# --- CONFIGURATION ---
pdb_folder = "/bucket/LaurinoU/PhD_student/Lilian/SequenceSimilarityNetwork/Sequences/ECODstructures/Structures/SAM_Binding_FaultyBAB/"
ecod_file = "/bucket/LaurinoU/PhD_student/Lilian/SequenceSimilarityNetwork/Sequences/ECODstructures/EcodDatabaseGrep/RestRossmannRelated/S_adenosyl_L_methionine.txt"
output_folder = "/bucket/LaurinoU/PhD_student/Lilian/SequenceSimilarityNetwork/Sequences/ECODstructures/Structures/Geo_Analysis/SAM_BABstructures_SecondRound2/"
os.makedirs(output_folder, exist_ok=True)

MIN_HELIX_LEN = 4
MIN_STRAND_LEN = 3
MIN_SHEET_SIZE = 4
MAX_BAB_LEN = 40  


# Helper function to extract numeric residue number from PDB ATOM line
def parse_resnum_from_line(line):
    resnum_str = line[22:26].strip()
    if not resnum_str: return None
    try: return int(resnum_str)
    except ValueError: return None


# --- Main loop ---
for pdbfile in os.listdir(pdb_folder):
    if not pdbfile.endswith(".pdb"): continue


    # Step 1: get pdb code
    basename = os.path.splitext(pdbfile)[0]  # e.g. "pdb2e5v"
    current_pdb = basename.replace("pdb", "").lower()  # -> "2e5v"
    print(f"\nProcessing {current_pdb}")

    
    
    # Step 2: find matching ECOD lines for the current pdb code in the file with all ECOD lines for domains within that superfamily (e.g. FAD/NAD(P) Rossmann fold binding proteins)
    rossmann_entries = []
    with open(ecod_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12: continue
            if parts[4].strip().lower() != current_pdb: continue
            chain_id = parts[5]
            domain_str = parts[6].split(",")[0]
            try:
                _, rng = domain_str.split(":")
                m = re.match(r'^\s*([+-]?\d+)[A-Za-z]*\s*-\s*([+-]?\d+)[A-Za-z]*\s*$', rng)
                if m:
                    rossmann_entries.append((chain_id, int(m.group(1)), int(m.group(2))))
            except: continue

    if not rossmann_entries: continue
    target_chain, ross_start, _ = sorted(rossmann_entries, key=lambda x: x[1])[0]

    
    
    #Step 3: Read PDB and save helix records, sheet records and all atom (and hetatm) lines
    helix_records, sheet_records, atom_lines = [], [], []
    with open(os.path.join(pdb_folder, pdbfile)) as f:
        for line in f:
            if line.startswith("HELIX"):
                c, s, e = line[19], int(line[21:25].strip()), int(line[33:37].strip())
                if (e - s + 1) >= MIN_HELIX_LEN: helix_records.append((c, s, e))
            elif line.startswith("SHEET"):
                sid, c, s, e = line[11:14].strip(), line[21], int(line[22:26].strip()), int(line[33:37].strip())
                if (e - s + 1) >= MIN_STRAND_LEN: sheet_records.append((c, sid, s, e))
            elif line.startswith(("ATOM", "HETATM")):
                new_line = "ATOM  " + line[6:] if line.startswith("HETATM") else line
                if line[17:20] == "MSE": new_line = new_line[:17] + "MET" + new_line[20:]
                atom_lines.append(new_line)

    
    
    #Step 4: Prepare search pool of possible first beta strands ordered by distance from the ECOD assigned start residue
    sheet_counts = Counter([sid for c, sid, s, e in sheet_records if c == target_chain])
    chain_helices = [(s, e) for c, s, e in helix_records if c == target_chain]
    
    # Sort all valid strands by distance from ECOD start
    candidate_strands = [
        (sid, s, e) for c, sid, s, e in sheet_records 
        if c == target_chain and sheet_counts[sid] >= MIN_SHEET_SIZE
    ]
    candidate_strands.sort(key=lambda x: abs(x[1] - ross_start))

    #Step 5: Find the most likely first beta strand through an iterative search to find the closest beta-alpha-beta segment to the ECOD assigned start that is shorter than 40 amino acids
    found_valid_bab = False
    
    for strand1 in candidate_strands:
        sid1, s1_start, s1_end = strand1
        
        # 1. Find the first helix after this specific strand1
        helices_after = [h for h in chain_helices if h[0] > s1_end]
        if not helices_after: continue
        h_start, h_end = min(helices_after, key=lambda x: x[0])
        
        # 2. Find the first strand after this helix
        all_strands = [(s, e) for c, sid, s, e in sheet_records if c == target_chain]
        strands_after = [s for s in all_strands if s[0] > h_end]
        if not strands_after: continue
        s2_start, s2_end = min(strands_after, key=lambda x: x[0])
        
        # 3. Check Length Constraint
        total_len = s2_end - s1_start
        if total_len <= MAX_BAB_LEN:
            
            #Step 6: Save the most likely beta-alpha-beta segment for this PDB
            region_start, region_end = max(1, s1_start - 3), s2_end + 3
            out_path = os.path.join(output_folder, f"BAB_{current_pdb}.pdb")
            with open(out_path, "w") as out:
                out.write(f"HELIX    1     1 {target_chain}{h_start:4d}    {target_chain:>1}{h_end:4d} 1\n")
                out.write(f"SHEET    1   A 2 {target_chain}{s1_start:4d}    {target_chain}{s1_end:4d}  0\n")
                out.write(f"SHEET    2   A 2 {target_chain}{s2_start:4d}    {target_chain}{s2_end:4d} -1\n")
                for line in atom_lines:
                    if line[21] == target_chain:
                        res = parse_resnum_from_line(line)
                        if res and region_start <= res <= region_end:
                            out.write(line)
            found_valid_bab = True
            break # Exit the search for this PDB file

    if not found_valid_bab:
        print(f"  Failed: No BAB segment under {MAX_BAB_LEN} aa found for {current_pdb}")
