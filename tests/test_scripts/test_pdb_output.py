#!/usr/bin/env python3

"""
Assess the completeness of a structure by comparing it to a reference structure.
Assumes that csymmatch has been used first to account for any origin shifts.
Outputs a JSON file with total, built and sequenced protein and nucleic acid residues.
Protein completeness is assessed using N, CA and C.
Nucleic completeness is assessed using C1', C2', C3', C4' and O4'.
For a reference residue to be classed as built,
all the atoms must be within a 1A radius of the same atom in the built structure.
For the residue to be classed as sequenced it must also be the same type.
"""

import argparse
import functools
import json
from dataclasses import dataclass
import dataclasses
import gemmi


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure")
    parser.add_argument("reference")
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("expected")
    return parser.parse_args()


def _main():
    args = _parse_args()
    structure = gemmi.read_structure(args.structure)
    reference = gemmi.read_structure(args.reference)

    search = gemmi.NeighborSearch(structure, max_radius=args.radius)
    search.populate(include_h=False)

    matching = functools.partial(_matching_residue, structure, search, args.radius)

    stats = _calculate_stats(reference, matching)

    if stats.nucleic_built >= int(args.expected):
        exit(0)
    else:
        exit(-1);


@dataclass
class Stats:
    protein_total: int = 0
    protein_built: int = 0
    protein_sequenced: int = 0
    nucleic_total: int = 0
    nucleic_built: int = 0
    nucleic_sequenced: int = 0


def _calculate_stats(reference, matching):
    stats = Stats()

    for chain in reference[0]:
        for residue in chain.first_conformer():
            info = gemmi.find_tabulated_residue(residue.name)
            if info.kind == gemmi.ResidueInfoKind.AA:
                atoms = {"N", "CA", "C"}
                if any(name not in residue for name in atoms):
                    continue

                stats.protein_total += 1

                if matching(residue, atoms, same_type=False):
                    stats.protein_built += 1

                if matching(residue, atoms, same_type=True):
                    stats.protein_sequenced += 1

            elif info.kind in (gemmi.ResidueInfoKind.RNA, gemmi.ResidueInfoKind.DNA):
                atoms = {"C1'", "C2'", "C3'", "C4'", "O4'"}

                if any(name not in residue for name in atoms):
                    continue
                stats.nucleic_total += 1

                if matching(residue, atoms, same_type=False):
                    stats.nucleic_built += 1

                if residue.name != "U" or "O2" in residue:  # UNK is U in Nautilus
                    if matching(residue, atoms, same_type=True):
                        stats.nucleic_sequenced += 1
    return stats


def _matching_residue(structure, search, radius, residue, atoms, same_type):
    return all(
        _matching_atom(structure, search, radius, residue, name, same_type)
        for name in atoms
    )


def _matching_atom(structure, search, radius, residue, name, same_type):
    for atom_alt in residue[name]:
        for mark in search.find_atoms(atom_alt.pos, "\0", radius=radius):
            cra = mark.to_cra(structure[0])
            if cra.atom.name == name:
                if not same_type or cra.residue.name == residue.name:
                    # print(gemmi.calculate_current_rmsd())
                    return True
    return False


if __name__ == "__main__":
    _main()
