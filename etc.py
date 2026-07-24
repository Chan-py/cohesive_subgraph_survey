import csv

def dict_to_csv(dict, key_col, value_col, file_name):
    with open(file_name, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([key_col, value_col])
        for k, v in dict.items():
            writer.writerow([k, v])
    print(f"File saved to {file_name}")