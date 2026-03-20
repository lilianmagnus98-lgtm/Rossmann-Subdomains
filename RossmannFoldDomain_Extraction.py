import os
import re

PDB_FOLDER = "/Path/to/PDB/Folder/"
ECOD_FILE = "/Path/to/ecod.latest.domain-DownloadFromECOD.txt"
OUTPUT_FOLDER = "/Path/to/Output/Folder/"


os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ---------------------------------------------------------------
# Parse ECOD and extract Rossmann fold region per PDB
# ---------------------------------------------------------------

def parse_ecod(ecod_file):
    """Return a dict: pdb -> list of (chain, [(start,end), (start,end)]) of each Rossmann fold listed for the pdb, including multiple chain"""
    raw_domains = {}    # pdb -> list of (chain, [(start,end)])

    with open(ecod_file) as f:
        for line in f:
            
            if "FAD/NAD(P)-binding domain" not in line:   #Define the Rossmann fold domain in question
                continue
            
            cols = line.strip().split("\t")
            if len(cols) < 7:
                continue

            pdb = cols[4]      # PDB code (e.g. "2e5v")
            chain = cols[5]    # Chain letter (e.g. "A")
            ranges = cols[6]   # e.g. "A:1-205,A:316-376"

            # Extract ranges like A:1-205,A:316-376
            region_list = []
            for r in ranges.split(","):
                m = re.match(r"([A-Za-z0-9]):(\d+)-(\d+)", r)
                if m:
                    ch, start, end = m.group(1), int(m.group(2)), int(m.group(3))
                    if ch == chain:     # sanity check
                        region_list.append((start, end))
            
            if not region_list:
                continue
            
            raw_domains.setdefault(pdb, []).append((chain, region_list))
            
    # -----------------------------------------------
    # Filter to keep only:
    #   • one chain per PDB (lowest alphabetically)
    #   • all domains from that chain
    # -----------------------------------------------
    final_domains = {}

    for pdb, entries in raw_domains.items():
        # collect chains that appear
        chains = sorted(list(set([c for (c, _) in entries])))

        # keep only the *first* chain alphabetically
        chosen_chain = chains[0]

        # filter entries for that chain only
        filtered = [(c, regions) for (c, regions) in entries if c == chosen_chain]

        final_domains[pdb] = filtered  # multiple entries if there are two Rossmann fold domains on one chain

    return final_domains

rossmann_map = parse_ecod(ECOD_FILE)
print(f"Found {len(rossmann_map)} Rossmann-containing PDB entries.")


# ---------------------------------------------------------------
# Helper: check if residue number lies in the extracted ranges
# ---------------------------------------------------------------

def residue_in_ranges(resnum, ranges):
    for start, end in ranges:
        if start <= resnum <= end:
            return True
    return False


# ---------------------------------------------------------------
# Process each PDB file and copy all sheet records, helix records and atom (and hetatm) lines within the specified range
# ---------------------------------------------------------------

for pdbcode, domain_list in rossmann_map.items():

    pdb_file = os.path.join(PDB_FOLDER, f"pdb{pdbcode}.pdb")
    if not os.path.isfile(pdb_file):
        continue

    # Multiple domains for the same PDB on the same chain
    for idx, (chain, ranges) in enumerate(domain_list, start=1):

        output_file = os.path.join(
            OUTPUT_FOLDER,
            f"{pdbcode}_{idx}_FADRoss.pdb"
        )

        with open(pdb_file) as fin, open(output_file, "w") as fout:
            for line in fin:
                record = line[:6].strip()

                # ATOM / HETATM
                if record in ("ATOM", "HETATM"):
                    line_chain = line[21].strip()
                    resnum = int(line[22:26])

                    if line_chain == chain and residue_in_ranges(resnum, ranges):
                        fout.write(line)

                # HELIX
                elif record == "HELIX":
                    hchain = line[19].strip()
                    hstart = int(line[21:25])
                    hend = int(line[33:37])
                    
                    if hchain == chain and (
                        residue_in_ranges(hstart, ranges)
                        or residue_in_ranges(hend, ranges)
                    ):
                        fout.write(line)

                # SHEET
                elif record == "SHEET":
                    schain = line[21].strip()
                    sstart = int(line[22:26])
                    send = int(line[33:37])

                    if schain == chain and (
                        residue_in_ranges(sstart, ranges)
                        or residue_in_ranges(send, ranges)
                    ):
                        fout.write(line)

        print(f"✓ Extracted {pdbcode} domain {idx} on chain {chain}")


