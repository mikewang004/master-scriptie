#!/usr/bin/env python3
"""
Sort LAMMPS dump file by atom ID and reassign mol_ids based on polymer length
The mol column is the 2nd column (index 1)
Usage: python sort_dump.py input.dump output.dump polymer_length
"""

import sys

def sort_and_reassign_mol_ids(input_file, output_file, polymer_length):
    """
    Read a LAMMPS dump file, sort atoms by ID, and reassign mol_ids
    
    Args:
        input_file: path to input dump file
        output_file: path to output sorted dump file
        polymer_length: number of atoms per polymer chain
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Find the ITEM: ATOMS line
    atom_line_idx = None
    header_line = None
    for i, line in enumerate(lines):
        if line.startswith('ITEM: ATOMS'):
            atom_line_idx = i
            header_line = line.strip()
            break
    
    if atom_line_idx is None:
        raise ValueError("Could not find 'ITEM: ATOMS' in file")
    
    # The mol column is always the 2nd column in the data (index 1)
    # Format: id mol xu yu zu
    mol_idx = 1
    
    # Header lines (everything before ITEM: ATOMS inclusive)
    header = lines[:atom_line_idx + 1]
    
    # Data lines (atoms)
    data_lines = lines[atom_line_idx + 1:]
    
    # Parse atoms
    atoms = []
    for line in data_lines:
        if line.strip():  # skip empty lines
            parts = line.strip().split()
            if len(parts) >= 2:  # Need at least id and mol
                atom_id = int(parts[0])
                atoms.append((atom_id, line, parts))
    
    # Sort by atom ID
    atoms.sort(key=lambda x: x[0])
    
    # Reassign mol_ids based on polymer length
    # mol_id 1: atoms 1 to polymer_length
    # mol_id 2: atoms polymer_length+1 to 2*polymer_length
    # etc.
    reassigned_atoms = []
    for idx, (atom_id, original_line, parts) in enumerate(atoms, start=1):
        # Calculate new mol_id (1-indexed)
        new_mol_id = (idx - 1) // polymer_length + 1
        
        # Replace the mol_id (2nd column, index 1)
        parts[mol_idx] = str(new_mol_id)
        
        # Reconstruct the line with proper spacing
        new_line = ' '.join(parts) + '\n'
        reassigned_atoms.append((atom_id, new_line))
    
    # Write output
    with open(output_file, 'w') as f:
        # Write header
        f.writelines(header)
        
        # Write reassigned atoms in sorted order
        for _, line in reassigned_atoms:
            f.write(line)
    
    num_molecules = (len(atoms) + polymer_length - 1) // polymer_length
    print(f"Sorted and reassigned {len(atoms)} atoms by ID")
    print(f"Polymer length: {polymer_length} atoms per chain")
    print(f"Number of molecules: {num_molecules}")
    print(f"Output written to: {output_file}")
    
    # Show first few and last few assignments as verification
    print("\nFirst 10 atom assignments:")
    for i in range(min(10, len(reassigned_atoms))):
        atom_id, line = reassigned_atoms[i]
        parts = line.strip().split()
        print(f"  Atom ID: {atom_id:6d} -> mol_id: {parts[mol_idx]}")
    
    if len(reassigned_atoms) > 10:
        print("  ...")
        print(f"Last 5 atom assignments:")
        for i in range(max(0, len(reassigned_atoms)-5), len(reassigned_atoms)):
            atom_id, line = reassigned_atoms[i]
            parts = line.strip().split()
            print(f"  Atom ID: {atom_id:6d} -> mol_id: {parts[mol_idx]}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python sort_dump.py <input_file> <output_file> <polymer_length>")
        print("Example: python sort_dump.py input.dump sorted.dump 100")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        polymer_length = int(sys.argv[3])
        if polymer_length <= 0:
            raise ValueError("Polymer length must be positive")
    except ValueError as e:
        print(f"Error: Invalid polymer length. {e}")
        sys.exit(1)
    
    try:
        sort_and_reassign_mol_ids(input_file, output_file, polymer_length)
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()