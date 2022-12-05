import matplotlib.pyplot as plt
import os
import json

def view_0_vs_7_completeness_graph():

    test_completeness_path = f'tests/test_angle_step/completeness_files/'

    all_data = []

    for file in os.scandir(test_completeness_path):
        with open(file.path, 'r') as json_file:
            json_data = json.load(json_file)

            data_tmp = {
                "file_name": int(file.name.split("_")[2].split('.')[0]),
                "nucleic_total": json_data["nucleic_total"],
                "nucleic_built": json_data["nucleic_built"],
                "nucleic_sequenced": json_data["nucleic_sequenced"]
            }

            all_data.append(data_tmp)

    all_data = sorted(all_data, key=lambda d: d['file_name'])

    names = [x['file_name'] for x in all_data]
    total = [x['nucleic_total'] for x in all_data]
    built = [(100*x['nucleic_built'])/x['nucleic_total'] for x in all_data]
    sequenced = [x['nucleic_sequenced'] for x in all_data]

    completeness_0 = [(100*x/y) for x, y in zip(built, total)]

    # plt.plot(names, total, label="Total Nucleic Acids")
    plt.plot(names, built, label="Nucleic Acids Built")
    # plt.plot(names,sequenced, label="Nucleic Acids Sequenced")
    plt.xlabel("Angle Step / degrees")
    plt.ylabel("Model complteness / %")
    plt.legend()
    plt.savefig('./tests/test_angle_step/output.png')

if __name__ == "__main__":
    view_0_vs_7_completeness_graph()