from Bio import SeqIO

def find_double_peak_sequences(fasta_file):
    double_peak_sequences = []
    double_peak_nucleotides = {'Y', 'y', 'M', 'm', 'K', 'k','S', 's', 'W', 'w', 'R', 'r'}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        if any(nucleotide in sequence for nucleotide in double_peak_nucleotides):
            double_peak_sequences.append(record)
    
    return double_peak_sequences

def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

# Directory of the fasta files
fasta_file = "../_data/_processed_data/473RT_mafft.fas"
output_file = "../_data/_processed_data/UniqueRT_473_doublepeak_.fas"
double_peak_sequences = find_double_peak_sequences(fasta_file)

# Write the double peak sequences to the output file
write_sequences_to_fasta(double_peak_sequences, output_file)

# Print the double peak sequences to the console
for seq_record in double_peak_sequences:
    print(f">{seq_record.id}")
    print(seq_record.seq)