import multiprocessing
import os
import subprocess
from dataclasses import dataclass

import gemmi
from tqdm import tqdm


@dataclass
class Parameters:
    test_list_dir: str = './tests/test_structures'
    data_dir: str = './tests/test_search_function'
    number_of_cycles: str = "3"
    library_list_path: str = './tests/test_files/rebuilt_filenames.txt'
    library_dir_path: str = './tests/test_library/'


def get_test_files():
    pdb_list = []

    for file in os.scandir(params.test_list_dir):
        if 'pdb' in file.name:
            prefix = file.name.split('.')[0]
            pdb_code = prefix.split('pdb')[1]
            if pdb_code not in pdb_list:
                pdb_list.append(pdb_code)
        else:
            pdb_code = file.name.split('.')[0]
            if pdb_code not in pdb_list:
                pdb_list.append(pdb_code)

    return pdb_list


def run_nautilus(test_pdb):
    out_file_path = os.path.join(params.data_dir, "pdb_files", f"{test_pdb}.pdb")
    test_pdb_mtz = os.path.join(params.test_list_dir, f"{test_pdb}.mtz")
    test_pdb_seq = os.path.join(params.test_list_dir, f"{test_pdb}.seq")

    library_export = "export LD_LIBRARY_PATH=$CLIB:$LD_LIBRARY_PATH"
    nautilus_cmd = f"""./cnautilus \
            -seqin {test_pdb_seq} \
            -mtzin {test_pdb_mtz} \
            -colin-fo FP,SIGFP \
            -colin-hl sfcalc.ABCD.A,sfcalc.ABCD.B,sfcalc.ABCD.C,sfcalc.ABCD.D \
            -cycles {params.number_of_cycles} \
            -anisotropy-correction \
            -pdblistin {params.library_list_path} \
            -pdblistdir {params.library_dir_path} \
            -title {test_pdb} \
            -pdbout {out_file_path}"""
    # -colin-free FREE \

    os.system(library_export)
    subprocess.run(nautilus_cmd, shell=True, stdout=subprocess.DEVNULL)


def worker(test_pdb):
    # run_nautilus(test_pdb)

    reference_structure_path = os.path.join(params.test_list_dir, f"pdb{test_pdb}.ent")
    search_function_structure_path = f"debug/sugar_positions/after_find_function_{test_pdb}.pdb"
    try:
        run_eval(reference_structure_path, search_function_structure_path, test_pdb)
    except:
        print("Error with ", test_pdb)

def main():
    test_list = get_test_files()

    with multiprocessing.Pool() as pool:
        x = list(tqdm(pool.imap_unordered(worker, test_list), total=len(test_list)))


def run_eval(path_to_ref, path_to_search_pdb, test_pdb):

    reference = gemmi.read_structure(path_to_ref)
    built = gemmi.read_structure(path_to_search_pdb)

    ref_neighbour_search = gemmi.NeighborSearch(reference, max_radius=5)
    ref_neighbour_search.populate(include_h=False)

    total = 0
    found = 0
    out_structure = gemmi.Structure()

    for chain in built[0]:
        for index, residue in enumerate(chain.first_conformer()):
            total += 1

            ref_P = residue.find_atom("P", "*")
            ref_O5 = residue.find_atom("O5'", "*")

            marks_p = ref_neighbour_search.find_atoms(ref_P.pos, radius=5)
            marks_o5 = ref_neighbour_search.find_atoms(ref_O5.pos, radius=5)

            if len(marks_p) > 0 and len(marks_o5) > 0:
                model = gemmi.Model(str(index + 1))

            alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

            for alpha_index, (p, o5) in enumerate(zip(marks_p, marks_o5)):
                cra_p = p.to_cra(reference[0])
                cra_o = o5.to_cra(reference[0])

                chain = gemmi.Chain(alphabet[alpha_index])

                if cra_p.atom.name == "P" and cra_o.atom.name == "O5'":
                    bond_length = cra_p.atom.pos.dist(cra_o.atom.pos)

                    vec_p = gemmi.Position(ref_P.pos.x - ref_O5.pos.x, ref_P.pos.y - ref_O5.pos.y,
                                           ref_P.pos.z - ref_O5.pos.z)
                    vec_o = gemmi.Position(cra_p.atom.pos.x - cra_o.atom.pos.x, cra_p.atom.pos.y - cra_o.atom.pos.y,
                                           cra_p.atom.pos.z - cra_o.atom.pos.z)

                    angle = vec_p.dot(vec_o)
                    angle_def = angle * (180 / 3.142)

                    if 1.65 > bond_length > 1.55 and angle_def < 20:
                        cra_p_res = cra_p.residue
                        chain.add_residue(residue, 1)
                        found += 1

                model.add_chain(chain)
                out_structure.add_model(model)

    out_structure.write_minimal_pdb(f'./debug/sugar_positions/gemmi_search_return_{test_pdb}.pdb')
    print(f"{test_pdb} - Built {found} out of {total} NAs ({(100 * found / total):2f} %)")


if __name__ == "__main__":
    params = Parameters()
    main()
