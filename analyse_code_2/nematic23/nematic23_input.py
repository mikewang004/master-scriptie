#!/usr/bin/env python3
# Usage: python convert_lammps.py <input_file> <Lchain> <output_file>

import sys
import os
import subprocess
import shutil
from pathlib import Path


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


def run_file(filename: str, workdir: str = "."):
    """
    Run ./nematic23 [filename] in workdir, move any newly created files
    to /output, and delete the input file.

    Returns a list of Paths to the moved files in /output.
    """
    workdir = Path(workdir).resolve()
    input_path = (workdir / filename).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Ensure /output exists
    output_dir = workdir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Snapshot files before running nematic23
    before_files = {p for p in workdir.iterdir() if p.is_file()}

    # Run ./nematic23 [filename] in the given working directory
    try:
        subprocess.run(
            ["./nematic23", input_path.name],
            cwd=workdir,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        # Process failed; propagate the error
        raise RuntimeError(f"nematic23 failed with exit code {e.returncode}") from e

    # Snapshot files after running nematic23
    after_files = {p for p in workdir.iterdir() if p.is_file()}

    # Newly created files (exclude the original input file and the executable)
    new_files = [
        p for p in (after_files - before_files)
        if p != input_path and p.name != "nematic23"
    ]

    moved_files = []
    for src in new_files:
        dst = output_dir / src.name
        # If a file with same name exists in /output, overwrite it
        if dst.exists():
            dst.unlink()
        shutil.move(str(src), str(dst))
        moved_files.append(dst)

    # Delete the original input file
    try:
        input_path.unlink()
    except FileNotFoundError:
        pass  # already deleted or moved elsewhere

    return moved_files

def main():
    # Lchain = 100
    # file_name = "equil_t_088_tdot_e-3_time_0.txt"
    # path_to_file = "../../../data/pva-100/quick_quench/%s" %(file_name)
    Lchain = 1000
    file_name = "PVA-%i_equil_t_088_tdot_e-3_run2_time_36000000.txt" %Lchain
    path_to_file = "../../../data/PVA-%i/equil/%s" %(Lchain, file_name)
    convert(path_to_file, Lchain, "./%s" %(file_name))
    run_file(file_name)


if __name__ == "__main__":
    main()


