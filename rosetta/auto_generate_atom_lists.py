#!/usr/bin/env python3
"""
Auto-generate decoy and linker atom alignment lists using Maximum Common Substructure (MCS).

This script automatically determines which atoms in the linker correspond to atoms
in the warheads (CBN and POI) by finding the Maximum Common Substructure between
the linker molecule and each warhead molecule. It eliminates the need for manual
specification of atom names, making the pipeline work with any linker PDB file.

Usage:
    # Basic usage - generate both decoy and linker atom lists
    python auto_generate_atom_lists.py \\
        --linker linkers/DY-0417_1.pdb \\
        -cbn CBN.pdb \\
        -poi POI.pdb \\
        -o auto_decoy_atom_list.txt \\
        -o_linker auto_linker_atom_list.txt

    # Also generate the delete atom list (all atoms in CBN and POI)
    python auto_generate_atom_lists.py \\
        --linker linkers/DY-0417_1.pdb \\
        -cbn CBN.pdb \\
        -poi POI.pdb \\
        -o auto_decoy_atom_list.txt \\
        -o_linker auto_linker_atom_list.txt \\
        --delete auto_decoy_atom_list_delete.txt

Dependencies:
    - RDKit (with MCS support: rdkit.Chem.rdFMCS)
"""

import argparse
import os
import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS


def load_molecule(pdb_file, remove_hs=True, mol2_file=None):
    """Load a PDB file and return an RDKit molecule with proper bond information.

    Uses proximity bonding to infer bonds, which is necessary for PDB files
    that lack CONECT records. For warhead molecules, a mol2 file can be
    provided to get correct bond information (recommended for SMARTS matching).

    When a mol2 file is provided, the function reads both files:
    - The mol2 file provides the correct chemical structure (bonds, formal charges)
    - The PDB file provides the atom names
    - Atoms are matched between the two files by 3D coordinates

    Args:
        pdb_file: Path to PDB file
        remove_hs: Whether to remove hydrogen atoms (recommended for MCS matching)
        mol2_file: Optional mol2 file with correct bond information

    Returns:
        RDKit Mol object (with PDB atom names from the PDB file, but
        bond information from the mol2 file if provided), or None if loading fails
    """
    if mol2_file and os.path.exists(mol2_file):
        # Read mol2 for correct bond information
        mol_mol2 = Chem.MolFromMol2File(mol2_file, sanitize=True, removeHs=False)
        if mol_mol2 is None:
            print(f"  Warning: Could not read mol2 file {mol2_file}, falling back to PDB")
        else:
            # Read PDB for atom names
            mol_pdb = Chem.MolFromPDBFile(
                pdb_file, sanitize=False, removeHs=False, flavor=0,
                proximityBonding=False
            )
            if mol_pdb is None:
                print(f"  Warning: Could not read PDB file {pdb_file}, using mol2 only")
                mol = mol_mol2
            else:
                # Transfer PDB atom names onto mol2 structure by matching coordinates
                mol = transfer_pdb_names(mol_mol2, mol_pdb)
                if mol is None:
                    print(f"  Warning: Could not transfer PDB names, using mol2 only")
                    mol = mol_mol2

            if remove_hs:
                mol = Chem.RemoveHs(mol)
            return mol

    # Fallback: read PDB with proximity bonding
    mol = Chem.MolFromPDBFile(
        pdb_file,
        sanitize=True,
        removeHs=False,
        flavor=0,
        proximityBonding=True
    )
    if mol is None:
        print(f"Error: Could not read {pdb_file}", file=sys.stderr)
        return None

    if remove_hs:
        mol = Chem.RemoveHs(mol)
    return mol


def transfer_pdb_names(mol_mol2, mol_pdb, tolerance=0.5):
    """Transfer PDB atom names onto a mol2 structure by matching 3D coordinates.

    Both molecules must have the same number of heavy atoms in approximately
    the same positions. Atoms are matched by finding the closest PDB atom
    to each mol2 atom within the given tolerance.

    Args:
        mol_mol2: RDKit molecule from mol2 file (correct chemistry)
        mol_pdb: RDKit molecule from PDB file (correct atom names)
        tolerance: Maximum distance (Angstrom) for coordinate matching

    Returns:
        New RDKit molecule with mol2 chemistry and PDB atom names, or None
    """
    from rdkit import Geometry

    # Get conformers
    conf_mol2 = mol_mol2.GetConformer(0)
    conf_pdb = mol_pdb.GetConformer(0)

    # Build a mapping: for each mol2 atom, find the closest PDB atom
    # Use a KD-tree for efficient matching
    pdb_positions = []
    for i in range(mol_pdb.GetNumAtoms()):
        pos = conf_pdb.GetAtomPosition(i)
        pdb_positions.append(pos)

    # Create new molecule with mol2 structure but PDB names
    mol = Chem.RWMol(mol_mol2)

    for i in range(mol_mol2.GetNumAtoms()):
        pos_mol2 = conf_mol2.GetAtomPosition(i)

        # Find closest PDB atom
        min_dist = float('inf')
        closest_pdb = -1
        for j, pos_pdb in enumerate(pdb_positions):
            dist = pos_mol2.Distance(pos_pdb)
            if dist < min_dist:
                min_dist = dist
                closest_pdb = j

        if min_dist < tolerance and closest_pdb >= 0:
            # Transfer PDB info
            pdb_atom = mol_pdb.GetAtomWithIdx(closest_pdb)
            pdb_info = pdb_atom.GetPDBResidueInfo()
            if pdb_info:
                new_pdb = Chem.AtomPDBResidueInfo()
                new_pdb.SetName(pdb_info.GetName())
                new_pdb.SetResidueName(pdb_info.GetResidueName())
                new_pdb.SetResidueNumber(pdb_info.GetResidueNumber())
                new_pdb.SetChainId(pdb_info.GetChainId())
                new_pdb.SetIsHeteroAtom(pdb_info.GetIsHeteroAtom())
                new_pdb.SetOccupancy(pdb_info.GetOccupancy())
                new_pdb.SetTempFactor(pdb_info.GetTempFactor())
                mol.GetAtomWithIdx(i).SetMonomerInfo(new_pdb)

    # Update conformer positions from mol2 (they should be the same)
    result = mol.GetMol()
    # Ensure the conformer has the correct positions
    conf = result.GetConformer(0)
    for i in range(result.GetNumAtoms()):
        conf.SetAtomPosition(i, conf_mol2.GetAtomPosition(i))

    return result


def get_pdb_atom_info(mol, atom_idx):
    """Get PDB atom name and residue name for a given atom index.

    Args:
        mol: RDKit Mol object
        atom_idx: Index of atom in the molecule

    Returns:
        Tuple of (atom_name, residue_name) or (None, None) if not found
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    pdb_info = atom.GetPDBResidueInfo()
    if pdb_info is not None:
        return pdb_info.GetName().strip(), pdb_info.GetResidueName().strip()
    return None, None


def get_mcs_pairs(mol_a, mol_b, atom_names_a, atom_names_b,
                  res_name_a, res_name_b, timeout=30, seed_smarts=None):
    """Find Maximum Common Substructure between two molecules and return atom pairs.

    Optionally uses a SMARTS pattern to constrain the MCS to include specific
    chemical substructures. This is useful for guiding the MCS to find the
    correct binding region of a PROTAC linker.

    Args:
        mol_a: First RDKit molecule
        mol_b: Second RDKit molecule
        atom_names_a: List of PDB atom names for mol_a atoms
        atom_names_b: List of PDB atom names for mol_b atoms
        res_name_a: PDB residue name for mol_a
        res_name_b: PDB residue name for mol_b
        timeout: Timeout in seconds for MCS search
        seed_smarts: Optional SMARTS pattern that the MCS must contain.
                     If provided, the MCS result will be a superset of this
                     pattern. If the pattern is not found in both molecules,
                     a warning is printed and unconstrained MCS is used.

    Returns:
        List of tuples: [(res_name_a atom_name_a, res_name_b atom_name_b), ...]
        Empty list if no MCS found
    """
    mcs_kwargs = dict(
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=timeout,
    )

    # Handle seed SMARTS: validate and add to MCS parameters
    if seed_smarts:
        smarts_mol = Chem.MolFromSmarts(seed_smarts)
        if smarts_mol is None:
            print(f"  Warning: Could not parse seed SMARTS: {seed_smarts}")
        else:
            has_a = mol_a.HasSubstructMatch(smarts_mol)
            has_b = mol_b.HasSubstructMatch(smarts_mol)
            if has_a and has_b:
                mcs_kwargs['seedSmarts'] = seed_smarts
                print(f"    Using seed SMARTS constraint: {seed_smarts}")
            else:
                missing = []
                if not has_a: missing.append(res_name_a)
                if not has_b: missing.append(res_name_b)
                print(f"  Warning: seed SMARTS not found in {'/'.join(missing)}, "
                      f"falling back to unconstrained MCS")

    mcs_result = rdFMCS.FindMCS([mol_a, mol_b], **mcs_kwargs)

    if mcs_result.numAtoms == 0:
        print(f"  Warning: No common substructure found between {res_name_a} and {res_name_b}")
        return []

    mcs_smarts = mcs_result.smartsString
    print(f"  MCS: {mcs_result.numAtoms} atoms, {mcs_result.numBonds} bonds")

    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    if mcs_mol is None:
        print(f"  Warning: Could not parse MCS SMARTS: {mcs_smarts}")
        return []

    match_a = mol_a.GetSubstructMatch(mcs_mol)
    match_b = mol_b.GetSubstructMatch(mcs_mol)

    if len(match_a) != len(match_b):
        print(f"  Warning: Match length mismatch: {len(match_a)} vs {len(match_b)}")
        return []

    pairs = []
    for idx_a, idx_b in zip(match_a, match_b):
        name_a = atom_names_a[idx_a] if idx_a < len(atom_names_a) else None
        name_b = atom_names_b[idx_b] if idx_b < len(atom_names_b) else None
        if name_a and name_b:
            pairs.append((f"{res_name_a} {name_a}", f"{res_name_b} {name_b}"))
        else:
            print(f"  Warning: Could not get atom names for indices {idx_a} -> {idx_b}")

    return pairs


def get_atom_names_from_mol(mol):
    """Extract PDB atom names from an RDKit molecule.

    Args:
        mol: RDKit Mol object

    Returns:
        List of atom name strings, or empty strings if PDB info not available
    """
    names = []
    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is not None:
            names.append(pdb_info.GetName().strip())
        else:
            names.append("")
    return names


def get_residue_name_from_mol(mol):
    """Extract PDB residue name from the first atom of an RDKit molecule.

    Args:
        mol: RDKit Mol object

    Returns:
        Residue name string, or 'UNK' if not available
    """
    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is not None:
            return pdb_info.GetResidueName().strip()
    return "UNK"


def group_by_rings_and_chains(mol, atom_indices):
    """Group atoms by ring systems and chain connectivity.

    This creates chemically meaningful groups:
    - Atoms in the same ring system are grouped together
    - Chain atoms are grouped by connectivity (3-4 atoms per group max)
    - Single atoms are kept as individual groups

    Args:
        mol: RDKit Mol object
        atom_indices: List of atom indices to group

    Returns:
        List of lists, each containing atom indices in a group
    """
    atom_set = set(atom_indices)
    if not atom_set:
        return []

    # Get ring information
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()

    # Find which rings overlap with our atom set
    # Process rings by size (largest first) to avoid atom duplication
    ring_groups = []
    ring_atoms_used = set()

    # Sort rings by size descending so larger rings are processed first
    sorted_rings = sorted(rings, key=len, reverse=True)

    for ring in sorted_rings:
        # Only include atoms that haven't been used yet
        ring_atoms = [a for a in ring if a in atom_set and a not in ring_atoms_used]
        if len(ring_atoms) >= 3:  # Only consider rings with at least 3 matched atoms
            ring_groups.append(ring_atoms)
            ring_atoms_used.update(ring_atoms)

    # For non-ring atoms, group by connectivity with max size limit
    non_ring_atoms = [a for a in atom_indices if a not in ring_atoms_used]
    non_ring_set = set(non_ring_atoms)

    # Group non-ring atoms by connectivity, but limit group size
    chain_groups = []
    visited = set()

    for start_idx in non_ring_atoms:
        if start_idx in visited:
            continue

        # BFS to find connected non-ring atoms
        group = []
        queue = [start_idx]
        visited.add(start_idx)

        while queue and len(group) < 4:  # Limit group size
            current = queue.pop(0)
            group.append(current)

            if len(group) >= 4:
                break

            atom = mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx in non_ring_set and n_idx not in visited:
                    visited.add(n_idx)
                    queue.append(n_idx)

        if group:
            chain_groups.append(sorted(group))

    # Handle any remaining unvisited non-ring atoms
    remaining = [a for a in non_ring_atoms if a not in visited]
    for rem in remaining:
        chain_groups.append([rem])

    # Combine ring groups and chain groups
    all_groups = ring_groups + chain_groups

    # Filter out empty groups
    all_groups = [g for g in all_groups if g]

    # Make sure every atom is in exactly one group
    grouped_atoms = set()
    for g in all_groups:
        grouped_atoms.update(g)
    ungrouped = [a for a in atom_indices if a not in grouped_atoms]
    for a in ungrouped:
        all_groups.append([a])

    return all_groups


def generate_atom_lists(linker_mol, warhead_mol, warhead_name,
                        linker_names, warhead_names,
                        linker_resname, warhead_resname,
                        seed_smarts=None):
    """Generate alignment atom pairs between a linker and a warhead molecule.

    Uses MCS to find common substructure, then groups matching atoms
    by connectivity. Optionally accepts a SMARTS pattern to constrain
    the MCS to include specific chemical substructures.

    Args:
        linker_mol: RDKit Mol of the linker
        warhead_mol: RDKit Mol of the warhead (CBN or POI)
        warhead_name: Name of the warhead for display ('CBN' or 'POI')
        linker_names: PDB atom names of linker atoms
        warhead_names: PDB atom names of warhead atoms
        linker_resname: PDB residue name of the linker
        warhead_resname: PDB residue name of the warhead
        seed_smarts: Optional SMARTS pattern to constrain the MCS

    Returns:
        Tuple of (warhead_groups, linker_groups) where each group is a list of
        atom name strings, or (None, None) if no MCS found
    """
    print(f"  Finding MCS between linker and {warhead_name}...")

    pairs = get_mcs_pairs(
        linker_mol, warhead_mol,
        linker_names, warhead_names,
        linker_resname, warhead_resname,
        seed_smarts=seed_smarts
    )

    if not pairs:
        return None, None

    # Build mapping from warhead atom idx to linker atom idx
    warhead_atoms = []
    linker_atoms = []

    # pairs[0] = linker key (e.g., "UNK N3"), pairs[1] = warhead key (e.g., "CBN N4")
    for linker_key, warhead_key in pairs:
        lk_name = linker_key.split()[-1]    # PDB atom name from linker
        wh_name = warhead_key.split()[-1]   # PDB atom name from warhead

        for idx, name in enumerate(linker_names):
            if name == lk_name:
                linker_atoms.append((lk_name, idx))
                break
        for idx, name in enumerate(warhead_names):
            if name == wh_name:
                warhead_atoms.append((wh_name, idx))
                break

    # Create mapping: warhead_idx -> linker_idx
    idx_map = {}
    for (_, wh_idx), (_, lk_idx) in zip(warhead_atoms, linker_atoms):
        idx_map[wh_idx] = lk_idx

    # Group by ring systems and chain connectivity in the warhead molecule
    warhead_indices = [idx for _, idx in warhead_atoms]
    wh_groups = group_by_rings_and_chains(warhead_mol, warhead_indices)

    # Build output groups
    warhead_groups = []
    linker_groups = []

    for group in wh_groups:
        wh_group = [f"{warhead_resname} {warhead_names[idx]}" for idx in group]
        lk_group = [f"{linker_resname} {linker_names[idx_map[idx]]}" for idx in group]
        warhead_groups.append(wh_group)
        linker_groups.append(lk_group)

    return warhead_groups, linker_groups


def write_atom_list_file(filename, groups):
    """Write atom list file in the format required by ternary_model_prediction.py.

    Each line has the format: "RESNAME ATOM1 ATOM2 ATOM3 ..."
    Atoms within a line are grouped together (centroid will be computed).

    Note: No header/comments are written because the existing read_atom_file()
    parser does not support comment lines.

    Args:
        filename: Output file path
        groups: List of groups, each group is a list of "RESNAME ATOMNAME" strings
    """
    with open(filename, 'w') as f:
        # Note: No header comments - the existing read_atom_file() parser
        # does not support comment lines. Each line must be a data line.
        for group in groups:
            if not group:
                continue
            first = group[0]
            resname = first.split()[0]
            atom_names = [item.split()[-1] for item in group]
            for item in group[1:]:
                rn = item.split()[0]
                if rn != resname:
                    print(f"  Warning: Mixed residue names in group: {resname} vs {rn}")
            line = f"{resname} " + " ".join(atom_names)
            f.write(line + "\n")


def generate_delete_list(cbn_pdb_file, poi_pdb_file, output_file):
    """Generate a list of all atoms in CBN and POI that should be deleted.

    Reads the original PDB files (with hydrogens) to get all atom names.

    Args:
        cbn_pdb_file: Path to CBN.pdb
        poi_pdb_file: Path to POI.pdb
        output_file: Output file path
    """
    cbn_full = Chem.MolFromPDBFile(cbn_pdb_file, sanitize=False, removeHs=False,
                                   proximityBonding=False)
    poi_full = Chem.MolFromPDBFile(poi_pdb_file, sanitize=False, removeHs=False,
                                   proximityBonding=False)

    with open(output_file, 'w') as f:
        for mol in [cbn_full, poi_full]:
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                pdb_info = atom.GetPDBResidueInfo()
                if pdb_info:
                    resname = pdb_info.GetResidueName().strip()
                    atomname = pdb_info.GetName().strip()
                    f.write(f"{resname} {atomname}\n")


def main():
    parser = argparse.ArgumentParser(
        description="""Auto-generate decoy and linker atom alignment lists for PROTAC ternary
                       complex modeling using Maximum Common Substructure (MCS).

                       This script eliminates the need to manually specify atom names for
                       alignment. It automatically finds the common chemical substructures
                       between the linker and warhead molecules (CBN and POI) and generates
                       the atom list files in the format required by ternary_model_prediction.py."""
    )

    parser.add_argument("-l", "--linker", required=True,
                        help="Linker PDB file (representative conformer)")
    parser.add_argument("-cbn", "--cbn", default=None,
                        help="CBN warhead PDB file (e.g., CBN.pdb)")
    parser.add_argument("-poi", "--poi", default=None,
                        help="POI warhead PDB file (e.g., POI.pdb)")
    parser.add_argument("-o", "--output_decoy", default="auto_decoy_atom_list.txt",
                        help="Output file for decoy atom alignment list")
    parser.add_argument("-o_linker", "--output_linker", default="auto_linker_atom_list.txt",
                        help="Output file for linker atom alignment list")
    parser.add_argument("--delete", default=None,
                        help="Output file for atoms to delete after merging (optional)")
    parser.add_argument("--timeout", type=int, default=30,
                        help="Timeout in seconds for MCS search (default: 30)")
    parser.add_argument("--cbn_smarts", default=None,
                        help="SMARTS pattern that the CBN-linker MCS must contain. "
                             "Guides the MCS to find the correct CBN-binding region. "
                             "Example: '[#7]1-[#6]-[#6]-[#6](-[#6]-1=[#8])=[#8]' for glutarimide")
    parser.add_argument("--poi_smarts", default=None,
                        help="SMARTS pattern that the POI-linker MCS must contain. "
                             "Guides the MCS to find the correct POI-binding region. "
                             "Example: '[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]#[#7]' for benzonitrile")
    parser.add_argument("--cbn_mol2", default=None,
                        help="CBN mol2 file with correct bond information (optional). "
                             "When provided, this is used for MCS matching instead of the PDB file, "
                             "giving correct chemical structure. Atom names are taken from the PDB file.")
    parser.add_argument("--poi_mol2", default=None,
                        help="POI mol2 file with correct bond information (optional). "
                             "When provided, this is used for MCS matching instead of the PDB file, "
                             "giving correct chemical structure. Atom names are taken from the PDB file.")

    args = parser.parse_args()

    if not args.cbn or not args.poi:
        print("Error: Both -cbn and -poi files are required", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("Auto-Generating Atom Alignment Lists")
    print("=" * 60)

    print(f"\n[1] Loading linker: {args.linker}")
    linker_mol = load_molecule(args.linker, remove_hs=True)
    if linker_mol is None:
        sys.exit(1)
    linker_names = get_atom_names_from_mol(linker_mol)
    linker_resname = get_residue_name_from_mol(linker_mol)
    print(f"    Linker: {linker_mol.GetNumAtoms()} heavy atoms, residue: {linker_resname}")

    print(f"\n[2] Loading warheads:")
    print(f"    CBN: {args.cbn}")
    if args.cbn_mol2:
        print(f"         (mol2: {args.cbn_mol2})")
    print(f"    POI: {args.poi}")
    if args.poi_mol2:
        print(f"         (mol2: {args.poi_mol2})")
    cbn_mol = load_molecule(args.cbn, remove_hs=True, mol2_file=args.cbn_mol2)
    poi_mol = load_molecule(args.poi, remove_hs=True, mol2_file=args.poi_mol2)
    if cbn_mol is None or poi_mol is None:
        sys.exit(1)

    cbn_names = get_atom_names_from_mol(cbn_mol)
    poi_names = get_atom_names_from_mol(poi_mol)
    cbn_resname = get_residue_name_from_mol(cbn_mol)
    poi_resname = get_residue_name_from_mol(poi_mol)
    print(f"    CBN: {cbn_mol.GetNumAtoms()} heavy atoms, residue: {cbn_resname}")
    print(f"    POI: {poi_mol.GetNumAtoms()} heavy atoms, residue: {poi_resname}")

    print(f"\n[3] Finding common substructures via MCS...")

    cbn_groups, linker_cbn_groups = generate_atom_lists(
        linker_mol, cbn_mol, "CBN",
        linker_names, cbn_names,
        linker_resname, cbn_resname,
        seed_smarts=args.cbn_smarts
    )

    poi_groups, linker_poi_groups = generate_atom_lists(
        linker_mol, poi_mol, "POI",
        linker_names, poi_names,
        linker_resname, poi_resname,
        seed_smarts=args.poi_smarts
    )

    print(f"\n[4] Writing output files...")

    decoy_groups = []
    linker_groups = []

    if cbn_groups:
        print(f"    CBN-Linker alignment: {len(cbn_groups)} groups")
        decoy_groups.extend(cbn_groups)
        linker_groups.extend(linker_cbn_groups)

    if poi_groups:
        print(f"    POI-Linker alignment: {len(poi_groups)} groups")
        decoy_groups.extend(poi_groups)
        linker_groups.extend(linker_poi_groups)

    if not decoy_groups:
        print("Error: No alignment pairs found!", file=sys.stderr)
        sys.exit(1)

    write_atom_list_file(args.output_decoy, decoy_groups)
    print(f"    Decoy atom list: {args.output_decoy}")

    write_atom_list_file(args.output_linker, linker_groups)
    print(f"    Linker atom list: {args.output_linker}")

    if args.delete:
        print(f"\n[5] Generating delete atom list: {args.delete}")
        generate_delete_list(args.cbn, args.poi, args.delete)
        print(f"    Delete atom list: {args.delete}")

    print(f"\n{'=' * 60}")
    print(f"Done! Generated {len(decoy_groups)} alignment groups.")
    print(f"")
    print(f"To use these files with ternary_model_prediction.py:")
    print(f"  python ../ternary_model_prediction.py \\")
    print(f"    -la {args.output_linker} \\")
    print(f"    -da {args.output_decoy} \\")
    print(f"    -dl decoy_list.txt \\")
    print(f"    -ll linker_list.txt \\")
    if args.delete:
        print(f"    -wd {args.delete} \\")
    print(f"    -c 0.4 \\")
    print(f"    -t default \\")
    print(f"    -r rmsd.txt")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()