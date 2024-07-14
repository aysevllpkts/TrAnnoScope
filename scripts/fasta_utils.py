# fasta_utils.py
import os

def count_sequences(fasta_file):
    total_sequences = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                total_sequences += 1
    return total_sequences

def find_fasta_with_max_sequences(fasta_files):
    max_sequences = 0
    max_fasta_file = ""
    for fasta_file in fasta_files:
        num_sequences = count_sequences(fasta_file)
        if num_sequences > max_sequences:
            max_sequences = num_sequences
            max_fasta_file = fasta_file
    return max_fasta_file, max_sequences

def calculate_chunks(num_sequences, chunk_size=5000):
    return (num_sequences + chunk_size - 1) // chunk_size


"""
# Functions
def count_sequences(fasta_file):
    total_sequences = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                total_sequences += 1
    return total_sequences

def find_fasta_with_max_sequences(fasta_files):
    max_sequences = 0
    max_fasta_file = ""
    for fasta_file in fasta_files:
        num_sequences = count_sequences(fasta_file)
        if num_sequences > max_sequences:
            max_sequences = num_sequences
    return max_sequences

def calculate_chunks(num_sequences, chunk_size=5000):
    return (num_sequences + chunk_size - 1) // chunk_size

# List of FASTA files
fasta_files = [os.path.join(OUTDIR, f"clustered/{sample}.clustered.fasta") for sample in SAMPLES]

# Find the FASTA file with the highest number of sequences
max_sequences = find_fasta_with_max_sequences(fasta_files)

# Calculate the number of chunks needed
chunk_size = config["blastn"]["split_size"]
chunks = calculate_chunks(max_sequences, chunk_size)

# Print the results (for debugging purposes)
print(f"Number of chunks needed: {chunks}")

##Indices for the fasta splited files, in the form of '001', '002'...
indexes=[]
for x in range(1,chunks+1): 
    indexes.append('%03d' % x)
print(indexes)
"""