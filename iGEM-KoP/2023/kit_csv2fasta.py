import csv

with open("./igem23kit.csv", newline="") as f:
    reader = csv.open(f)
    for row in reader:
        print(row)
