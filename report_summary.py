import os
import pandas as pd

data_dir = "results/reports/"
output_file = "combined_report.tsv"

relevant_columns = ["vs2_score", "vs2_hallmark", "vs2_viral", "genomad_viral", "genomad_hallmarks"]

def process_reports(data_dir, output_file):
    combined_data = []

    for file in os.listdir(data_dir):
        if file.endswith("_report.tsv"):
            file_path = os.path.join(data_dir, file)
            sample_name = file.split("_")[0]  + '_' + file.split("_")[1] # change if your sample names contains '_'

            try:
                report = pd.read_csv(file_path, sep="\t")
            except Exception as e:
                print(f"Error reading {file}: {e}")
                continue

            filtered_rows = report[
                report[relevant_columns].notna().any(axis=1) |
                report["base_taxonomy"].str.contains("Viruses", na=False)
            ]

            filtered_rows.insert(0, "Sample_Name", sample_name)
            combined_data.append(filtered_rows)

    if combined_data:
        final_data = pd.concat(combined_data, ignore_index=True)
        final_data.to_csv(output_file, sep="\t", index=False)
        print(f"Combined filtered report saved to {output_file}")
    else:
        print("No valid data to combine.")

process_reports(data_dir, output_file)

