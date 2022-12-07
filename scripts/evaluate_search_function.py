import gemmi


def main():

    reference_structure_path = './data/1hr2/1hr2_final.pdb'
    search_function_structure_path = './debug/sugar_positions/after_find_function.pdb'

    reference = gemmi.read_structure(reference_structure_path)
    built = gemmi.read_structure(search_function_structure_path)

    ref_neighbour_search = gemmi.NeighborSearch(reference, max_radius=5)
    ref_neighbour_search.populate(include_h=False)
    
    total = 0
    found = 0
    out_structure = gemmi.Structure()

    for chain in built[0]:
        for index, residue in enumerate(chain.first_conformer()):
            total+=1
            
            ref_P = residue.sole_atom("P")
            ref_O5 = residue.sole_atom("O5'")
        
            marks_p = ref_neighbour_search.find_atoms(ref_P.pos, radius=5)
            marks_o5 = ref_neighbour_search.find_atoms(ref_O5.pos, radius=5)

            if len(marks_p)>0 and len(marks_o5)>0:
                model = gemmi.Model(str(index+1))

            alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuwxyz"

            for index, (p, o5) in enumerate(zip(marks_p, marks_o5)):
                cra_p = p.to_cra(reference[0])
                cra_o = o5.to_cra(reference[0])

                chain = gemmi.Chain(alphabet[index])

                if cra_p.atom.name == "P" and cra_o.atom.name == "O5'":
                    bond_length = cra_p.atom.pos.dist(cra_o.atom.pos)

                    vec_p = gemmi.Position(ref_P.pos.x-ref_O5.pos.x, ref_P.pos.y-ref_O5.pos.y, ref_P.pos.z-ref_O5.pos.z)
                    vec_o = gemmi.Position(cra_p.atom.pos.x-cra_o.atom.pos.x, cra_p.atom.pos.y-cra_o.atom.pos.y, cra_p.atom.pos.z-cra_o.atom.pos.z)

                    angle = vec_p.dot(vec_o)
                    angle_def = angle * (180/3.142)

                    if bond_length < 1.61 and bond_length > 1.58 and angle_def < 20 :
                        cra_p_res = cra_p.residue
                        chain.add_residue(residue, 1)
                        found+=1
               
                model.add_chain(chain)
                out_structure.add_model(model)

    out_structure.write_minimal_pdb('./debug/sugar_positions/gemmi_search_return.pdb')
    print(f"Built {found} out of {total} NAs ({(100*found/total):2f} %)")
    

    # print(x)
if __name__ == "__main__":
    main()