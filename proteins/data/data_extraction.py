import pandas as pd
import os

from typing import List, Dict

results_folder = "./proteins/results"
data_folder = "./proteins/sequences"

all_data = []


nucleocapsid_proteins = ["nucleocapsid protein", "nucleocapsid phosphoprotein", "N protein"]
spike_proteins = ["spike protein", "spike glycoprotein", "S protein"]

def extract_accessions(
        data_folder: str, 
        target_proteins: List[List[str]]
) -> Dict[str, Dict[str, List[str]]]:
    """
    Get all the accessions for `target_proteins` from the .fasta sequences in the `data_folder`.
    `target_proteins` is a list of synonims for a given protein, and we use the first element of
    that list to identify the protein.
    """
    all_accessions = {}

    for filename in os.listdir(data_folder):
        file_path = os.path.join(data_folder, filename)
        virus_name = filename.split(".")[0]
        all_accessions[virus_name] = {}

        with open(file_path, "r") as f:
            for line in f:
                if line[0] == ">":
                    i = line.find("|")
                    j = line.find("[")

                    for target_protein in target_proteins:
                        if line[(i + 1) : (j - 1)] in target_protein:
                            if target_protein[0] not in all_accessions[virus_name]:
                                all_accessions[virus_name][target_protein[0]] = []
                            all_accessions[virus_name][target_protein[0]].append(line[1 : (i - 1)])

    return all_accessions

def find_key(dict: Dict[str, List[str]], value) -> str:
    return next((k for k, v in dict.items() if value in v), None)
                    
def create_data(
        results_folder: str,
        all_accessions: Dict[str, Dict[str, List[str]]],
        output_file: str
):

    for filename in os.listdir(results_folder):
        file_path = os.path.join(results_folder, filename)

        parts = filename.replace(".out", "").split(".")
        repeat_type = parts[1].split("_")[0]
        virus_name = parts[0]

        if not os.path.isfile(file_path):
            continue

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()

                if line.count(",") == 7:
                    data = line.split(",")
                    protein = find_key(all_accessions[virus_name], data[0])
                    print

                    if protein:
                        data.extend([virus_name, protein, repeat_type])
                        all_data.append(data)

    df = pd.DataFrame(
        all_data, 
        columns = [
            "sequence_id",
            "left_start",
            "left_end",
            "right_start",
            "right_end",
            "length",
            "left_repeat",
            "right_repeat",
            "virus_name",
            "protein_name",
            "repeat_type"
        ]
    )

    df.to_csv(output_file, index = False)

if __name__ == "__main__":
    all_accessions = extract_accessions(data_folder = data_folder, target_proteins = [spike_proteins, nucleocapsid_proteins])
    create_data(results_folder = results_folder, all_accessions = all_accessions, output_file = "./proteins/data/virus_protein_repeats.csv")
