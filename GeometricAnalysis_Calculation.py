import os
import numpy as np


BAB_pdb_folder = "/path/to/pdb/folder/"
output_file = "/path/to/output/file.txt"



def safe_int_from_field(s):     #Helper function for parsing the correct start and end residues
    digits = ''.join(filter(str.isdigit, s))
    if not digits:
        raise ValueError(f"No digits in field: '{s}'")
    return int(digits)

def calculate_strand_vector_endpoints(atom_lines, start_res, end_res):  #finds the endpoints of any strand region and calculates the strand vector
    """
    For beta strands: calculate the centroid of the top two residues and the centroid of the bottom two residues for the C-terminal and N-terminal point, respectively. The strand vector goes through these two centroids.
    """
    ca_coords = []
    
    for line in atom_lines:
        if line.startswith(("ATOM","HETATM")):
            atom_name = line[12:16].strip()
            res_num = int(line[22:26].strip())

            if atom_name == "CA" and start_res <= res_num <= end_res:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ca_coords.append([x, y, z])


    ca_array = np.array(ca_coords)

    # Strands shorter than 4 aminoacids will be disregarded
    if len(ca_array) >= 4:
        strand_Nter = np.mean(ca_array[:2], axis=0) # Mean of first two residues is the N-terminal
        strand_Cter = np.mean(ca_array[-2:], axis=0) # Mean of last two residues is the C-terminal
        
        
        strand_vector = strand_Cter - strand_Nter
    else:
        return None, None, None 
    
    return strand_vector, strand_Nter, strand_Cter

def calculate_helix_vector_endpoints(atom_lines, start_res, end_res): #finds the central endpoints of any helix and calculates the helix vector
    """
    For helices:  calculate the centroid of the top three residues and the centroid of the bottom three residues for the N-terminal and C-terminal point, respectively. The helix vector goes through these two centroids.
    """
    ca_coords = []

    # Collecting all Ca coordinates in order
    for line in atom_lines:
        if line.startswith(("ATOM","HETATM")):
            atom_name = line[12:16].strip()
            res_num = int(line[22:26].strip())

            if atom_name == "CA" and start_res <= res_num <= end_res:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ca_coords.append([x, y, z])
    
    
    ca_array = np.array(ca_coords)

    # Helixes shorter than 5 aminoacids will be disregarded
    if len(ca_coords) >= 5:
        # Take bottom 3 and top 3 residues 
        Nter_three = ca_array[:3]   # first 3 residues is the N-terminal
        Cter_three = ca_array[-3:]  # last 3 residues is the C-terminal

        # Calculate center points (centroids) and vector
        helix_Nter_center = np.mean(Nter_three, axis=0)
        helix_Cter_center = np.mean(Cter_three, axis=0)

        helix_vector = helix_Nter_center - helix_Cter_center
    else:
        return None, None, None

    return helix_vector, helix_Nter_center, helix_Cter_center

def calculate_angle_between_vectors(vec1, vec2):
    """
    Calculate the angle between two vectors in degrees (0-90°). Used for S1H, S2H, and SHS angles.
    """
    if vec1 is None or vec2 is None:
        return None
    
    # Normalize vectors
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    
    if norm1 == 0 or norm2 == 0:
        return None
    
    vec1_norm = vec1 / norm1
    vec2_norm = vec2 / norm2
    
    # Calculate dot product 
    dot_product = np.dot(vec1_norm, vec2_norm)
    dot_product = np.clip(dot_product, -1.0, 1.0)
    
    # Calculate angle from the cosine angle/dot product
    angle_rad = np.arccos(dot_product)
    angle_deg = np.degrees(angle_rad)
    
    # Return acute angle (0-90 degrees)
    
    if angle_deg > 90:
        angle_deg = 180 - angle_deg
    
    return angle_deg
    

    
def endpoint_distances(strand_Nter, strand_Cter, helix_Nter_center, helix_Cter_center):
  
    '''
    Calculate distances between strand and helixes. Used for S1H Top Dist., S1H Bottom Dist., S2H Top Dist., and S2H Bottom Dist. 
    '''

    bottom_distance    = np.linalg.norm(strand_Nter - helix_Cter_center)
    top_distance = np.linalg.norm(strand_Cter - helix_Nter_center)

    if bottom_distance < top_distance:
        scenario = "V-shape"
    else:
        scenario =  "A-shape"

    return scenario, top_distance, bottom_distance



def SH_twist(strand_Nter, strand_Cter, helix_Nter_center, helix_Cter_center):
    
    '''
    Calculate the twist of one beta strand - alpha helix structure
    '''
    #vectors between consecutive points
    helix_vector = helix_Nter_center - helix_Cter_center
    strand_vector = strand_Cter - strand_Nter
    stoh_vector = strand_Cter - helix_Nter_center  #Vector between Strand C-terminal point and Helix N-terminal point

    #Normal vectors to the two planes
    plane1 = np.cross(strand_vector, stoh_vector)
    plane2 = np.cross(stoh_vector, helix_vector)

    #Check for degenerate cases (collinear points)
    if np.linalg.norm(plane1) == 0 or np.linalg.norm(plane2) == 0:
        return 0.0
    
    #Normalize the normal vectors and the central bond
    nplane1 = plane1 / np.linalg.norm(plane1)
    nplane2 = plane2 / np.linalg.norm(plane2)
    nstoh_vector = stoh_vector / np.linalg.norm(stoh_vector)

    #Calculate the twist/torsion angle. Gives angle from -180 to 180 degrees
    m1 = np.cross (nplane1, nstoh_vector)
    x = np.dot(nplane1, nplane2)
    y = np.dot(m1, nplane2)
    twist = np.degrees(np.arctan2(y, x))

    return twist


def SHS_Vectors(strand1_C, strand2_C, helix_Nter_center):
    
    '''
    Calculate vectors between C-terminal of strands  and N-terminal of the helix
    '''

    S1H_vector = helix_Nter_center - strand1_C
    S2H_vector = helix_Nter_center - strand2_C

    return S1H_vector, S2H_vector




# --- Main loop ---

'''
Output file:
S1_Sta: Start of strand 1
S1_End: End of strand 1
S2_Sta: Start of strand 2
S2_End: End of strand 2
H_Sta: Start of helix 
H_End: End of helix
S1H_Angle: Angle between strand 1 and the helix. Between 0 - 90 degrees
S2H_Angle: Angle between strand 2 and the helix. Between 0 - 90 degrees
S1H_VA: Whether the top or the bottom of strand 1 and the helix lean together. Can discriminate between "V" and "A" looking structures. (V = bottom closer, A = top closer)
S2H_VA: Whether the top or the bottom of strand 2 and the helix lean together
S1H_TD: Distance between the top centroid of strand 1 and the top centroid of the helix 
S1H_BD: Distance between the bottom centroid of strand 1 and the bottom centroid of the helix
S2H_TD: Distance between the top centroid of strand 2 and the top centroid of the helix
S2H_BD: Distance between the bottom centroid of strand 2 and the bottom centroid of the helix
GL_Len: Length of the glycine rich loop. From end of strand 1 to the start of the helix.
SHS_Angle: Angle betweeen the vector from the top of strand 1 to the top of the helix and the vector from the top of strand 2 and the top of the helix
SH1_Twi: Twist/Torsion Angle between strand 1 and the helix
SH2_Twi: Twist/Torsion Angle between strand 2 and the helix
'''

with open(output_file, "w") as out_f:
    out_f.write("PDB\tS1_Sta\tS1_End\tS2_Sta\tS2_End\tH_Sta\tH_End\tGL_Len\tS1H_Angle\tS2H_Angle\tS1H_VA\tS2H_VA\tS1H_TD\tS1H_BD\tS2H_TD\tS2H_BD\tSHS_Angle\tSH1_Twi\tSH2_Twi\n")
    for pdbfile in os.listdir(BAB_pdb_folder):
        if not pdbfile.endswith(".pdb"):
            continue

        # Step 1: get pdb code from the filename
        basename = os.path.splitext(pdbfile)[0]   # e.g. "pdb2e5v"
        current_pdb = basename.replace("BAB_", "").lower()  # -> "2e5v"
        print(f"\nProcessing {current_pdb}")

        #Step 2: define Strand and Helix regions

        helix_records, strand_records, atom_lines = [], [], []
        with open(os.path.join(BAB_pdb_folder, pdbfile)) as f:
            for line in f:
                if line.startswith("HELIX"):
                    helix_start = safe_int_from_field(line[20:26])    #start residue of one helix
                    helix_end = safe_int_from_field(line[32:38])       #end residue of one helix
                    helix_records.append((helix_start, helix_end))
            
                elif line.startswith("SHEET"):
                    strand_start = safe_int_from_field(line[20:27])  #start residue of one sheet
                    strand_end = safe_int_from_field(line[32:38])       #end residue of one sheet
                    strand_records.append((strand_start, strand_end))
                elif line.startswith(("ATOM","HETATM")):
                    atom_lines.append(line)

         #Skip files with too many or too few HELIX/SHEET records
        if len(helix_records) > 1:
            print(f"Skipping {pdbfile}: found {len(helix_records)} HELIX records (expected 1)")
            continue
            
        if len(strand_records) > 2:
            print(f"Skipping {pdbfile}: found {len(strand_records)} SHEET records (expected ≤2)")
            continue

        if len(strand_records) < 2:
            print(f"Not enough strands found in {pdbfile}")
            continue
    
        if len(helix_records) < 1:
            print(f"No helix found in {pdbfile}")
            continue
        
        #Step 3: Assign Strand1 and Strand2
        else:
            # sort by start residue and assign strands
            strand_records.sort(key=lambda x: x[0])

            strand1_start, strand1_end = strand_records[0]
            strand2_start, strand2_end = strand_records[1]


            #Calculate Glycine-Loop Lengths
            Gloop_Length = helix_start - strand1_end - 1
            
            
            
            # Step 4: Calculate vectors
            res_s1 = calculate_strand_vector_endpoints(atom_lines, strand1_start, strand1_end)
            res_s2 = calculate_strand_vector_endpoints(atom_lines, strand2_start, strand2_end)
            res_h  = calculate_helix_vector_endpoints(atom_lines, helix_start, helix_end)

            # Check if any returned None (as too short or had missing CA atoms)
            if res_s1[0] is None or res_s2[0] is None or res_h[0] is None:
                print(f"Skipping {current_pdb}: Secondary structure too short for reliable vector calculation.")
                continue

            # Unpack the returned values
            strand1_vector, strand1_Nter, strand1_Cter = res_s1
            strand2_vector, strand2_Nter, strand2_Cter = res_s2
            helix_vector, helix_Nter_center, helix_Cter_center = res_h
        
            
            
            
            # Step 5: Calculate S1H and S2H angles
            angle_S1H = calculate_angle_between_vectors(strand1_vector, helix_vector)
            angle_S2H = calculate_angle_between_vectors(strand2_vector, helix_vector)
        

            #Format angles for output (show "NA" if None)
            angle_S1H_str = f"{angle_S1H:.2f}" if angle_S1H is not None else "NA"
            angle_S2H_str = f"{angle_S2H:.2f}" if angle_S2H is not None else "NA"


            
            
            #Step 6: Calculate Endpoint Distances and V or A scenario
            strand1_scenario, strand1_top_dist, strand1_bottom_dist  = endpoint_distances(strand1_Nter, strand1_Cter, helix_Nter_center, helix_Cter_center)
            strand2_scenario, strand2_top_dist, strand2_bottom_dist  = endpoint_distances(strand2_Nter, strand2_Cter, helix_Nter_center, helix_Cter_center)
            
            #Format endpoint distances for output
            Strand1_top_distance = f"{strand1_top_dist:.2f}"
            Strand1_bottom_distance = f"{strand1_bottom_dist:.2f}"
            Strand2_top_distance = f"{strand2_top_dist:.2f}"
            Strand2_bottom_distance = f"{strand2_bottom_dist:.2f}"
            
            
            
            
            #Step 7: Calculate the twist/torsion angle of one beta strand alpha helix structure
            S1H_twist = SH_twist(strand1_Nter, strand1_Cter, helix_Nter_center, helix_Cter_center)
            S2H_twist = SH_twist(strand2_Nter, strand2_Cter, helix_Nter_center, helix_Cter_center)

            #Format twist/torsion angles for output
            S1H_Twist = f"{S1H_twist:.2f}"
            S2H_Twist = f"{S2H_twist:.2f}"
            
            
            
            
            #Step 8: Calculate the SHS angle 
            S1toH_vector, S2toH_vector = SHS_Vectors(strand1_Cter, strand2_Cter, helix_Nter_center)
            S1HS2_angle = calculate_angle_between_vectors(S1toH_vector, S2toH_vector)
            SHS_Angle = f"{S1HS2_angle:.2f}"

            #Step 9: Write everything in file
            out_f.write(f"{current_pdb}\t{strand1_start}\t{strand1_end}\t{strand2_start}\t{strand2_end}\t{helix_start}\t{helix_end}\t{Gloop_Length}\t{angle_S1H_str}\t{angle_S2H_str}\t{strand1_scenario}\t{strand2_scenario}\t{Strand1_top_distance}\t{Strand1_bottom_distance}\t{Strand2_top_distance}\t{Strand2_bottom_distance}\t{SHS_Angle}\t{S1H_Twist}\t{S2H_Twist}\n")
            
            


