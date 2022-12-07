import gemmi


def main():

    reference_structure_path = './data/1hr2/1hr2_final.pdb'
    search_function_structure_path = './debug/sugar_positions/after_find_function.pdb'

    reference = gemmi.read_structure(reference_structure_path)
    built = gemmi.read_structure(search_function_structure_path)

    ref_neighbour_search = gemmi.NeighborSearch(reference, max_radius=3)
    ref_neighbour_search.populate(include_h=False)

    pos = gemmi.Position(56,67,78)
    marks = ref_neighbour_search.find_atoms(pos, "\0", radius=3)
    # print(marks)
    # quit()
    # print(built)
    tar_neighbour_search = gemmi.NeighborSearch(built[0], built.cell, 3)
    # print(built[0])
    total = 0
    found = 0
    for chain in built[0]:
        for residue in chain.first_conformer():
            total+=1
            ref_P = residue.sole_atom("P")

            marks = ref_neighbour_search.find_atoms(ref_P.pos, radius=3)
            for mark in marks:
                cra = mark.to_cra(reference[0])
                if cra.atom.name == "P":
                    # print(cra.atom)
                    found+=1

    print(f"{found=}, {total=}")
    print(f"{100*found/total}%")

    # print(x)
if __name__ == "__main__":
    main()