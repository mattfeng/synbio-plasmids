from Bio import SeqIO

# Path to your GenBank file
genbank_file = "example.gbk"

# Read the GenBank file
for record in SeqIO.parse(genbank_file, "genbank"):
    print(f"Record ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print(f"Description: {record.description}")
    print(f"Length: {len(record.seq)}")

enzymes = {
    "BsaI": (6, 1, 5), # restriction site len,
    "BsmBI": (),

}
indices = {
    "A": (75, 918),
    "B": (67, 910)
}

overhangs = {
    "1": "GGAG",
    "2": "TACT",
    "3": "AATG",
    "4": "GCTT",
    "5": "CGCT",
}