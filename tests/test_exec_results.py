import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import json


def view_results():

    test_completeness_path = f'tests/test_comparisons/base_nautilus/completeness_files'
    alternate_dir = 'tests/test_output/correlation_function/1hr2_library/completeness_files'

    all_data = []
    for file in os.scandir(test_completeness_path):

        file_name = file.name.split('.')[0]

        file_data = {
            "file_name": file_name
        }
        with open(file.path, 'r') as json_file:
            json_data = json.load(json_file)

            self_test_data = { 
                "nucleic_total": json_data["nucleic_total"],
                "nucleic_built": json_data["nucleic_built"],
                "nucleic_sequenced": json_data["nucleic_sequenced"]
            }   

            file_data["self_test"] = self_test_data

        alternate_path = os.path.join(alternate_dir, f"{file_name}.json")
        if os.path.isfile(alternate_path):
            with open(alternate_path, 'r') as json_file:
                json_data = json.load(json_file)

                sole_test_data = { 
                    "nucleic_total": json_data["nucleic_total"],
                    "nucleic_built": json_data["nucleic_built"],
                    "nucleic_sequenced": json_data["nucleic_sequenced"]
                }

                file_data["sole_test"] = sole_test_data
        else: 
            print(f"{file_name} not found!")

        if os.path.isfile(alternate_path):
            all_data.append(file_data)

    
    all_data = sorted(all_data,key=lambda x : x['self_test']['nucleic_total'], reverse=True)

    names = [x['file_name'] for x in all_data]
    total = [x["self_test"]['nucleic_total'] for x in all_data]
    self_built_0 = [x["self_test"]['nucleic_built'] for x in all_data]
    sole_built_3 = [x["sole_test"]['nucleic_built'] for x in all_data]

    completeness_0 = [(100*x/y) for x, y in zip(self_built_0, total)]
    completeness_3 = [(100*x/y) for x, y in zip(sole_built_3, total)]

    colors = []
    labels = []
    
    for name, x,y in zip(names, completeness_0, completeness_3):
        if x<y: 
            colors.append('green')
            labels.append('')
        elif x==y: 
            colors.append('blue')
            labels.append('')

        else: 
            colors.append('red')
            labels.append(name)

    delta = []

    for x,y in zip(sole_built_3, self_built_0):
        delta.append(x-y)
   
    gs = gridspec.GridSpec(2, 2)

    fig = plt.figure(figsize=(12,6))
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,:])

    ax0.plot(names, total, label="Total Nucleic Acids")
    ax0.plot(names, self_built_0, label="NA Built - Original Scoring")
    ax0.plot(names, sole_built_3, label="NA Built - New Scoring")
    ax1.bar(names, delta, label="NA Built - Delta", color=colors)
    # plt.plot(names, sole_built, label="NA Built - 1hr2 only")

    ax0.set_xlabel("PDB Test Structure")
    ax1.set_title("Difference in number of NA built")
    ax1.set_xlabel("PDB Test Structure")

    ax0.set_xticklabels(names, rotation=90)
    ax1.set_xticklabels(names, rotation=90)

    # plt.xticks(rotation=90)
    ax0.set_ylabel("Number")
    ax1.set_ylabel("Delta")

    ax0.legend()
    # ax1.legend() 

    ax2.scatter(completeness_0, completeness_3, c=colors)
    ax2.set_xlabel("Model completeness - Original Scoring / %")
    ax2.set_ylabel("Model completeness - New Scoring / %")
    
    for i, text in enumerate(labels): 
        ax2.annotate(text, (completeness_0[i], completeness_3[i]-10))

    lims = [
    np.min([ax2.get_xlim(), ax2.get_ylim()]),  # min of both axes
    np.max([ax2.get_xlim(), ax2.get_ylim()]),  # max of both axes
]
    ax2.plot(lims, lims, 'k-', alpha=0.75, zorder=0)

    plt.title("Orignal Scoring function vs New Scoring Function - Lib 1hr2 only")
    plt.tight_layout()
    plt.savefig(f'./tests/test_output/correlation_function/1hr2_library/updated_find+grow.png', dpi = 500)
    # plt.show()
if __name__ == "__main__":
    view_results()