from Bio import Entrez

# Set your email here
Entrez.email = "aysevilpektas@gmail.com"

# Function to fetch species name for a given taxid
def fetch_species_name(taxid):
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    return records[0]["ScientificName"]

# Path to the file with UniProt IDs and taxids
input_file = '../taxonomy_info/uniprot_taxid_results.txt'

# Read the file and fetch species names
with open(input_file, 'r') as f:
    results = []
    for line in f:
        uniprot_id, taxid = line.strip().split('\t')
        try:
            species_name = fetch_species_name(taxid)
            results.append((uniprot_id, taxid, species_name))
        except Exception as e:
            print(f"Error fetching species name for taxid {taxid}: {e}")
            results.append((uniprot_id, taxid, "N/A"))

# Save the results to a new file
output_file = '../taxonomy_info/uniprot_taxid_species_results.txt'
with open(output_file, 'w') as f:
    for uniprot_id, taxid, species_name in results:
        f.write(f"{uniprot_id}\t{taxid}\t{species_name}\n")

# Print the results
for uniprot_id, taxid, species_name in results:
    print(f"UniProt ID: {uniprot_id}, TaxID: {taxid}, Species Name: {species_name}")

