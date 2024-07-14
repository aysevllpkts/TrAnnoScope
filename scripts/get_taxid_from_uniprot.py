import gzip

# Function to extract UniProt ID and taxid
def extract_uniprot_taxid(file_path):
    with gzip.open(file_path, 'rt') as file:
        uniprot_id = None
        taxid = None
        results = []

        for line in file:
            if line.startswith('ID'):
                uniprot_id = line.split()[1]
            elif line.startswith('OX'):
                taxid = line.split('=')[1].split(';')[0]
            elif line.startswith('//'):
                if uniprot_id and taxid:
                    results.append((uniprot_id, taxid))
                uniprot_id = None
                taxid = None
        
        return results

# Path to your uniprot_sprot.dat.gz file
file_path = "/home/aysevil/MolGen/faststorage/aysevil/greenland_shark/FLTransAnnot/data/TRINOTATE_DATA_DIR/uniprot_sprot.dat.gz"
results = extract_uniprot_taxid(file_path)

# Print the results
for uniprot_id, taxid in results:
    print(f"UniProt ID: {uniprot_id}, TaxID: {taxid}")

# Save the results to a file
output_file = 'uniprot_taxid_results.txt'
with open(output_file, 'w') as f:
    for uniprot_id, taxid in results:
        f.write(f"{uniprot_id}\t{taxid}\n")

