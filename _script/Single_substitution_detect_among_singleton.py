from Bio import SeqIO
import csv

def read_fasta_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[str(record.seq)] = record.id
    return sequences

def find_one_mismatch_sequences(sequences1, sequences2):
    one_mismatch_sequences = []
    for seq1, id1 in sequences1.items():
        for seq2, id2 in sequences2.items():
            if len(seq1) == len(seq2) and sum(1 for a, b in zip(seq1, seq2) if a != b) == 1:
                one_mismatch_sequences.append((id1, id2))
    return one_mismatch_sequences

def write_to_csv(data, output_file):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Sequence1_ID", "Sequence2_ID"])
        writer.writerows(data)

# Directory of the fasta files
fasta_file_176 = '../_data/_processed_data/singleton176_mafft.fas'
fasta_file_473 = '../_data/_processed_data/unique473_mafft.fas'
output_csv = "../_data/_processed_data/Single_subsititution176.csv"

sequences_176 = read_fasta_sequences(fasta_file_176)
sequences_473 = read_fasta_sequences(fasta_file_473)

one_mismatch_sequences = find_one_mismatch_sequences(sequences_176, sequences_473)

write_to_csv(one_mismatch_sequences, output_csv)

print(f"Found {len(one_mismatch_sequences)} sequences with one mismatch.")

import csv

input_csv = "../_data/_processed_data/Single_subsititution176.csv"
unique_ids_csv = "../_data/_processed_data/Unique_sequence_IDs.csv"

sequence1_ids = set()
with open(input_csv, mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        sequence1_ids.add(row["Sequence1_ID"])  # keep only unique IDs in set

# Write the unique sequence IDs to the CSV file
with open(unique_ids_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Unique_Sequence_ID"])
    for uid in sorted(sequence1_ids):
        writer.writerow([uid])

# Print the number of unique sequence IDs to the console
print(f"Unique Sequence1_IDs (no duplicates): {len(sequence1_ids)}")