import argparse
from os.path import isfile, join
from os import listdir
import csv


def get_column_score(cols, mapping, col_name):
    if col_name in mapping:
        score = 0
        col = cols[mapping[col_name]]
        peps = col.split("|")
        for pep in peps:
            if "[" in pep and not "{d}" in pep:
                s = pep.split("[")[1]
                s = float(s.split("]")[0])
                score += s
        return score
    return -1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input", help="input files or directories", required=True, nargs="+")
    parser.add_argument("-mode", "--mode", help="'DIR' or 'FILE', depending on whether the downloaded PIRCHE job resulted in a single file or a director with several files", required=True)
    parser.add_argument("-output", "--output", help="input CSV", required=True)
    parser.add_argument("-donorIdentifier", "--donorIdentifier", help="string that is used in donor IDs to distinguish from patients (e.g. 'd' or 'donor')", required=True)

    args = parser.parse_args()

    files = []
    if args.mode == "DIR":
        for path in args.input:
            files.extend([ join(path,f) for f in listdir(path) if isfile(join(path,f))])
    if args.mode == "FILE":
        for path in args.input:
            files.append(path)
    print("Considering a total of " + str(len(files)) + " file(s)")

    with open(args.output, "w") as output:
        writer = csv.writer(output, delimiter=';')
        writer.writerow(["id", "A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "DRB1_1", "DRB1_2", "DRB3_1", "DRB3_2", "DRB4_1", "DRB4_2", "DRB5_1", "DRB5_2", "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2", "PIRCHE_II", "PIRCHE_II_corr", "PIRCHE_II_A", "PIRCHE_II_B", "PIRCHE_II_C", "PIRCHE_II_DRB1", "PIRCHE_II_DRB345", "PIRCHE_II_DQA1", "PIRCHE_II_DQB1", "PIRCHE_II_DPA1", "PIRCHE_II_DPB1", "PIRCHE_DRB1_A", "PIRCHE_DRB1_B", "PIRCHE_DRB1_C", "PIRCHE_DRB1_DRB1", "PIRCHE_DRB1_DRB345", "PIRCHE_DRB1_DQA1", "PIRCHE_DRB1_DQB1", "PIRCHE_DRB1_DPA1", "PIRCHE_DRB1_DPB1", "PIRCHE_DRB345_A", "PIRCHE_DRB345_B", "PIRCHE_DRB345_C", "PIRCHE_DRB345_DRB1", "PIRCHE_DRB345_DRB345", "PIRCHE_DRB345_DQA1", "PIRCHE_DRB345_DQB1", "PIRCHE_DRB345_DPA1", "PIRCHE_DRB345_DPB1", "PIRCHE_DQ_A", "PIRCHE_DQ_B", "PIRCHE_DQ_C", "PIRCHE_DQ_DRB1", "PIRCHE_DQ_DRB345", "PIRCHE_DQ_DQA1", "PIRCHE_DQ_DQB1", "PIRCHE_DQ_DPA1", "PIRCHE_DQ_DPB1", "PIRCHE_DP_A", "PIRCHE_DP_B", "PIRCHE_DP_C", "PIRCHE_DP_DRB1", "PIRCHE_DP_DRB345", "PIRCHE_DP_DQA1", "PIRCHE_DP_DQB1", "PIRCHE_DP_DPA1", "PIRCHE_DP_DPB1"])
        for f in files:
            data = {}
            mapping = {}
            with open(f, "r") as input:
                print("Read " + f)
                header = 2
                try:
                    for row in input:
                        cols = row.strip().split(",")
                        if header == 1:
                            for idx, k in enumerate(cols):
                                mapping[k] = idx
                        if not header > 0:
                            id = cols[0]
                            cols[0] = cols[0].replace("donor_", "")
                            if args.donorIdentifier in id:
                                data[id] = cols
                        header -= 1
                except UnicodeDecodeError:
                    print("Skip")



            for key in data:
                new_row = data[key][0:23]
                new_row.append(data[key][28])
                new_row.append(data[key][267])

                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_A_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_B_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_C_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB1_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB1_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB5_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DQA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DQB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DPA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DPB1_Epitopes"))

                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_A_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_B_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_C_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DRB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB1_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB1_Presents_DRB5_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DQA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DQB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DPA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DPB1_Epitopes"))

                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_A_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_B_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_C_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB5_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB5_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DQA1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DQA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DQB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DPA1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DPA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DRB3_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DPB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DPB1_Epitopes"))

                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_A_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_B_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_C_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB5_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DQA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DQB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DPA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DPB1_Epitopes"))

                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_A_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_B_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_C_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB3_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB4_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB5_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DQA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DQB1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DPA1_Epitopes"))
                new_row.append(get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DPB1_Epitopes"))

                writer.writerow(new_row)