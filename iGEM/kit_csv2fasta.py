import csv

with open("./igem23kit.csv", newline="") as fin:
    reader = csv.reader(fin)
    columns = {
        1: "part_id",
        3: "url",
        4: "part_type",
        16: "copy_number",
        17: "resistance"
    }
    for row in reader:
        data = {columns[i]: row[i] for i in columns}
        name = row[2]
        sequence = row[21]
        plate = row[8]
        well = row[9]

        header = ">" + "|".join(f"{k}={v}" for k, v in data.items())

        with open(f"./IG23KP_{plate}{well} {name.replace('/', '_')}.fasta", "w") as fout:
            fout.write(header + "\n")
            fout.write(sequence + "\n")
