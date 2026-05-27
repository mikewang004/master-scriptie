#!/usr/bin/env python3
# Usage: python convert_lammps.py <input_file> <Lchain> <output_file>

import sys
import os

def convert(input_file, Lchain, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Parse header
    n_atoms = None
    box = []
    atom_start = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "ITEM: NUMBER OF ATOMS":
            n_atoms = int(lines[i+1].strip())
            i += 2
        elif line == "ITEM: BOX BOUNDS pp pp pp":
            for j in range(3):
                box.append(lines[i+1+j].strip())
            i += 4
        elif line.startswith("ITEM: ATOMS"):
            atom_start = i + 1
            break
        else:
            i += 1

    if n_atoms is None or len(box) != 3 or atom_start is None:
        raise ValueError("Could not parse input file header.")

    # Derived quantities
    n_chains = n_atoms // Lchain
    n_bonds  = n_atoms - n_chains

    # Parse atoms: columns are id mol xu yu zu
    atoms = []
    for line in lines[atom_start:]:
        parts = line.split()
        if len(parts) < 5:
            continue
        orig_id = int(parts[0])
        xu, yu, zu = parts[2], parts[3], parts[4]
        atoms.append((orig_id, xu, yu, zu))

    # Sort by original id
    atoms.sort(key=lambda x: x[0])
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # Write output
    with open(output_file, 'w') as f:
        f.write("Nematic system data file\n")
        #f.write("\n")
        f.write(f"{n_atoms} atoms\n")
        f.write("1 atom types\n")
        f.write(f"{n_bonds} bonds\n")
        f.write("181 bond types\n")
        f.write("181 angles\n")
        f.write("1 angle types\n")
        #f.write("\n")
        f.write(f"{box[0]} xlo xhi\n")
        f.write(f"{box[1]} ylo yhi\n")
        f.write(f"{box[2]} zlo zhi\n")
        #f.write("\n")
        f.write("ITEM: ATOMS id mol xu yu zu\n")
        for new_id, (orig_id, xu, yu, zu) in enumerate(atoms, start=1):
            mol_id = (new_id - 1) // Lchain + 1
            f.write(f"{new_id} {mol_id} {xu} {yu} {zu}\n")

    print(f"Done. {n_atoms} atoms, {n_chains} chains, {n_bonds} bonds written to {output_file}")


def main():
    Lchain = 100
    file_name = "equil_t_088_tdot_e-3_time_0.txt"
    path_to_file = "../../../data/pva-100/quick_quench/%s" %(file_name)
    convert(path_to_file, Lchain, "input/%s" %(file_name))


if __name__ == "__main__":
    main()


