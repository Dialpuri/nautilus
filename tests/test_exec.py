import os
import csv
import shutil
from dataclasses import dataclass
import subprocess
from tqdm import tqdm
import multiprocessing


@dataclass
class Parameters:
    data_dir: str = './tests/test_output/correlation_function/1hr2_library/'
    test_list_dir: str = './tests/test_structures'
    number_of_cycles: str = "3"
    library_list_path: str = './tests/test_files/rebuilt_filenames.txt'
    library_dir_path: str = './tests/test_library/'


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
-pdbout {out_file_path}"""
# -colin-free FREE \
# -pdblistin {params.library_list_path} \
# -pdblistdir {params.library_dir_path} \

    os.system(library_export)
    subprocess.run(nautilus_cmd, shell=True, stdout=subprocess.DEVNULL)


def run_completeness_script(test_pdb):
    out_file_path = os.path.join(params.data_dir, "pdb_files", f"{test_pdb}.pdb")
    comparison_pdb_path = os.path.join(params.test_list_dir, f"pdb{test_pdb}.ent")
    completeness_out_path = os.path.join(params.data_dir, "completeness_files", f"{test_pdb}.json")
    python_cmd = f"python scripts/completeness.py {out_file_path} {comparison_pdb_path} {completeness_out_path}"

    if os.path.isfile(out_file_path):
        subprocess.run(python_cmd, shell=True, stdout=subprocess.DEVNULL)


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


def worker(test_pdb):
    run_nautilus(test_pdb)
    run_completeness_script(test_pdb)


def main():
    test_list = get_test_files()

    with multiprocessing.Pool() as pool:
        x = list(tqdm(pool.imap_unordered(worker, test_list), total=len(test_list)))

    # for test_pdb in tqdm(test_list):
    #     run_nautilus(test_pdb)
    #     run_completeness_script(test_pdb)
        # quit()

def generate_test_dir_structure():
    pdb_file_path = os.path.join(params.data_dir, "pdb_files")
    completeness_file_path = os.path.join(params.data_dir, "completeness_files")

    if os.path.isdir(pdb_file_path):
        shutil.rmtree(pdb_file_path)

    os.mkdir(pdb_file_path)

    if os.path.isdir(completeness_file_path):
        shutil.rmtree(completeness_file_path)

    os.mkdir(completeness_file_path)


if __name__ == "__main__":
    params = Parameters()
    generate_test_dir_structure()
    main()
