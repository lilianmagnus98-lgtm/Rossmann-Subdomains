import os
import re
from collections import defaultdict


# Helper function to ensure correct extraction of residue numbers and handling insertion codes
def parse_residue(res_str, icode=''):
    match = re.match(r'(\d+)([A-Za-z]?)', res_str)
    if not match:
        raise ValueError(f'Invalid residue number: {res_str}')
    num, ic = match.groups()
    return (int(num), ic or ' ')

#Helper function to sort residues correctly
def res_cmp_key(res):
    return (res[0], res[1])



def process_pdb(pdb_path, out_dir_first, out_dir_second):
    """
    Identifies a Rossmann-fold beta-sheet topology and splits the PDB into two halves.
    
    The function searches for a 3-strand Rossmann core (beta-alpha-beta-alpha-beta) 
    within a single chain. It identifies the 'start' strand and a 'third' parallel 
    strand to determine a structural breakpoint.

    Args:
        pdb_path: Path to the input .pdb file.
        out_dir_first: Directory to save the N-terminal half (up to the 3rd beta strand).
        out_dir_second: Directory to save the C-terminal half (from the 4th beta strand).

    Returns:
        bool: True if a valid Rossmann topology was found and files were written, 
              False otherwise.

    Raises:
        RuntimeError: If a SHEET line is malformed or cannot be parsed.
    
    Outputs:
        - Writes {base}_FirHalf.pdb to out_dir_first.
        - Writes {base}_SecHalf.pdb to out_dir_second.
    """


    with open(pdb_path) as f:
        lines = [line for line in f if line.startswith(('ATOM', 'HETATM', 'SHEET', 'HELIX'))]

    #Step 1: Parse SHEET lines to extract fixed-width columns
    sheets = []
    for line in lines:
        if line.startswith('SHEET'):
            try:
                sheet_id = line[11:15].strip()  # sheet label (e.g., AA2)
                strand_num = int(line[7:10].strip())
                strand_count = int(line[14:16].strip())
                chain = line[21].strip()
                init_res = line[22:27].strip() # first residue of the strand
                end_res = line[33:38].strip() # last residue of the strand
                sense = int(line[38:40].strip()) # sense vs antisense (0 or 1 in PDB files)
                sheets.append({
                    'line': line,
                    'sheet_id': sheet_id,
                    'strand_num': strand_num,
                    'strand_count': strand_count,
                    'chain': chain,
                    'init': parse_residue(init_res),
                    'end': parse_residue(end_res),
                    'sense': sense
                })
            except Exception as e:
                raise RuntimeError(f'Error parsing SHEET line in {pdb_path}: {line.strip()}') from e
    
    #Step 2: Group and renumber strands per sheet (start at 1, contiguous) to find the sheet with the most strands in the domain 
    sheet_groups = defaultdict(list)
    for s in sheets:
        sheet_groups[s['sheet_id']].append(s)

    largest_sheet = None
    max_strands = 0
    for sid, group in sheet_groups.items():
        # Sort by original strand_num
        group_sorted = sorted(group, key=lambda x: x['strand_num'])
        # Renumber sequentially from 1
        for i, s in enumerate(group_sorted, 1):
            s['strand_num'] = i
        # Finding the larest sheet (most strands)
        if len(group_sorted) > max_strands and len(group_sorted) > 3:
            max_strands = len(group_sorted)
            largest_sheet = {'id': sid, 'strands': group_sorted}
    
    if not largest_sheet:
        print(f' {os.path.basename(pdb_path)}: no sheet >3 strands')
        return False

    strands = largest_sheet['strands']
    sheet_id = largest_sheet['id']

    

    #Precompute strand index in sheet (1-indexed → 0-indexed)
    strand_to_idx = {s['strand_num']: i for i, s in enumerate(strands)}

    # Helper: get sheet index of a strand (0-based)
    def sheet_idx_of(strand):
        return strand_to_idx[strand['strand_num']]


    # Step 3: Find the orientation of the strands in relation to strand 1 (parallel vs antiparallel). The Rossmann fold should have 5 or 6 parallel strands. 
    # Precompute absolute orientation for each strand in sheet (1 = parallel to sheet start, -1 = antiparallel).
    # Strand 1: orientation = +1 (reference)
    # Strand i (i>1): orientation[i] = orientation[i-1] * (1 if sense_i == 1 else -1)
    orientation = {}
    for s in strands:
        n = s['strand_num']
        if n == 1:
            orientation[n] = 1
        else:
            # Use current strand's sense (relation to previous)
            orientation[n] = orientation[n-1] * (1 if s['sense'] == 1 else -1)
        s['abs_orient'] = orientation[n]
            

    # Sequence-ordered strands (by init residue, same chain as sheet majority? — but we'll filter by chain later)
    seq_ordered = sorted(strands, key=lambda s: res_cmp_key(s['init']))



    # Step 4: Iterative search to find the first three strands of the Rossmann fold
    candidate = None
    for i in range(min(3, len(seq_ordered))):  
        s = seq_ordered[i]
        idx = strand_to_idx[s['strand_num']]
        # Condition 1: first strand of the Rossmann fold should not be at the edge of the sheet
        if idx == 0 or idx == len(strands) - 1:
            continue
        # Condition 2: next strand in *sequence* (same chain) must be adjacent in the same sheet as strand 1
        next_in_seq = None
        for s2 in seq_ordered[i+1:]:
            if s2['chain'] == s['chain']:
                next_in_seq = s2
                break
        if not next_in_seq:
            continue
        next_idx = strand_to_idx[next_in_seq['strand_num']]
        if abs(next_idx - idx) != 1:
            continue
        candidate = s
        break

    if not candidate:
        print(f'{os.path.basename(pdb_path)}: no valid Rossmann start strand')
        return False

    # Step 5: Find 3rd Rossmann strand = 2nd *later* strand in sequence, same sheet, same chain, AND same absolute orientation to strands 1 and 2
    later_strands = [
        s for s in seq_ordered
        if (s['chain'] == candidate['chain']
            and res_cmp_key(s['init']) > res_cmp_key(candidate['init'])
            and s['sheet_id'] == sheet_id
            and s['abs_orient'] == candidate['abs_orient'])   
    ]
    
    #later_strands[0] = second beta strand of the Rossmann fold. later_strands[1] = third beta strand of the Rossmann fold etc.
    if len(later_strands) < 2:
        print(f'{os.path.basename(pdb_path)}: <2 later strands with same orientation for 3rd beta') #There should be atleast two more beta strand in the same sheet after strand 1
        return False

    third_beta = later_strands[1]  #Third beta strand of the rossmann fold
    end_first_half = third_beta['end'][0] + 2  #The end of the first half is defined as two residues after the end of the third beta strand of the Rossmann fold


    # Step 4: Find the fourth strand of the Rossmannd fold (strand after linker)
    
    # We used later_strands[1] for the third beta strand
    # So later_strands[2] is the forth beta strand
    if len(later_strands) < 3:
        print(f'[UNPARSEABLE] {os.path.basename(pdb_path)}: need ≥3 later strands for topology 4 (got {len(later_strands)})') #Less than 4 strands in the fold
        return False

    fourth_beta = later_strands[2]  # fourth beta strand of the Rossmann fold
    start_second_half = fourth_beta['init'][0] - 2  #The start of the second half is defined as two residues before the start of the forth beta strand
    # ===========================================================


    # Step 5: Extract lines from the pdb if they are within the cut off 
    atom_lines_fir = [] #Atom lines within the First Half Subdomain
    atom_lines_sec = [] #Atom lines within the Second Half Subdomain
    ss_lines_fir = [] #Secondary structure within the First Half Subdomain
    ss_lines_sec = [] #Secondary structure within the Second Half Subdomain
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain = line[21]
            res_str = line[22:27].strip()
            try:
                resnum, icode = parse_residue(res_str)
            except:
                continue
            if chain == candidate['chain'] and resnum <= end_first_half: #The first half includes any atom lines between the start of the Rossmann fold domain and the 2 residues after the end of the third strand
                atom_lines_fir.append(line)
            elif chain == candidate['chain'] and resnum >= start_second_half: #The second half includes any atom line between the 2 residues before the start of the forth strand and the end of the Rossmann fold domain
                atom_lines_sec.append(line)
        elif line.startswith(('SHEET', 'HELIX')):
            # For SHEET/HELIX: include if *any* residue in record ≤ cutoff
            chain = line[21]
            # SHEET: init 22-27, end 33-38; HELIX: init 22-27, end 33-37
            try:
                init_res = parse_residue(line[22:27].strip())
                end_res = parse_residue(line[33:38].strip())
            except:
                continue
            if chain == candidate['chain'] and init_res[0] <= end_first_half:
                ss_lines_fir.append(line)
            elif chain == candidate['chain'] and init_res[0] >= start_second_half:
                ss_lines_sec.append(line)

    # Step 6:  Write output of First and Second Half Subdomain in new files in their respective folders
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    out_path_first = os.path.join(out_dir_first, f'{base}_FirHalf.pdb')
    out_path_second = os.path.join(out_dir_second, f'{base}_SecHalf.pdb') 
    os.makedirs(out_dir_first, exist_ok=True)
    os.makedirs(out_dir_second, exist_ok=True)
    #First Half Subdomain
    with open(out_path_first, 'w') as f:
        f.writelines(ss_lines_fir)
        f.writelines(atom_lines_fir)
    #Second Half Subdomain
    with open(out_path_second, 'w') as f:
        f.writelines(ss_lines_sec)
        f.writelines(atom_lines_sec)

    return True

# ─────────────── MAIN ───────────────
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 4:
        print('Usage: python split_ros_by_beta.py <input_dir> <output_dir_1> <output_dir_2>')
        sys.exit(1)
    in_dir, out_dir_first, out_dir_second = sys.argv[1], sys.argv[2], sys.argv[3]
    for fname in os.listdir(in_dir):
        if fname.endswith('.pdb'):
            process_pdb(os.path.join(in_dir, fname), out_dir_first, out_dir_second)
