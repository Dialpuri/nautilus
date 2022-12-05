import os
import csv
from dataclasses import dataclass
import subprocess
from tqdm import tqdm
import multiprocessing
import numpy as np


@dataclass
class Parameters:
    data_dir: str = './tests/test_angle_step/'
    test_list_dir: str = './data/1hr2/'
    number_of_cycles: str = "3"


def run_nautilus(test_pdb, angle_step):
    out_file_path = os.path.join(params.data_dir, "pdb_files", f"{test_pdb}_{angle_step}.pdb")
    test_pdb_mtz = os.path.join(params.test_list_dir, f"{test_pdb}.mtz")
    test_pdb_seq = os.path.join(params.test_list_dir, f"{test_pdb}.seq")

    library_export = "export LD_LIBRARY_PATH=$CLIB:$LD_LIBRARY_PATH"
    nautilus_cmd = f"""./cnautilus \
-seqin {test_pdb_seq} \
-mtzin {test_pdb_mtz} \
-colin-fo FP,SIGFP \
-colin-fc FWT,PHWT \
-colin-free FREE \
-cycles {params.number_of_cycles} \
-anisotropy-correction \
-search-step {angle_step} \
-pdbout {out_file_path}"""

    os.system(library_export)
    subprocess.run(nautilus_cmd, shell=True, stdout=subprocess.DEVNULL)


def run_completeness_script(test_pdb, angle_step):
    out_file_path = os.path.join(params.data_dir, "pdb_files", f"{test_pdb}_{angle_step}.pdb")
    comparison_pdb_path = os.path.join(params.test_list_dir, f"{test_pdb}.pdb")
    completeness_out_path = os.path.join(params.data_dir, "completeness_files", f"{test_pdb}_{angle_step}.json")
    python_cmd = f"python scripts/completeness.py {out_file_path} {comparison_pdb_path} {completeness_out_path}"

    if os.path.isfile(out_file_path):
        subprocess.run(python_cmd, shell=True, stdout=subprocess.DEVNULL)


def worker(angle_step):
    test_pdb = '1hr2_final'
    run_nautilus(test_pdb, angle_step)
    run_completeness_script(test_pdb, angle_step)


def main():
    list_of_angles = np.arange(1, 20, 2)

    with multiprocessing.Pool() as pool:
        x = list(tqdm(pool.imap_unordered(worker, list_of_angles), total=len(list_of_angles)))

    # test_pdb = '1hr2_final'
    # for angle_step in tqdm(list_of_angles):
    #     run_nautilus(test_pdb, angle_step)
    #     run_completeness_script(test_pdb, angle_step)


if __name__ == "__main__":
    params = Parameters()

    main()
