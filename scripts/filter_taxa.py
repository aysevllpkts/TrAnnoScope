import sys
import json
import pandas as pd
from functools import reduce

def filter_identifiers(data, output_file, taxa_filter):
    """Filter identifiers based on taxa and save to a text file.

    Args:
        data (str): Path to the input data file (BlobTools table).
        output_file (str): Path to the output text file.
        taxa_filter (dict): Dictionary containing filter parameters for taxa removal.

    Returns:
        str: Message indicating the status of the operation.
    """
    try:
        # Load the BlobTools table
        data_df = pd.read_csv(data, sep='\t')

        # Extract the filter parameters for the 'phylum' or other groups
        phylum_filter = taxa_filter.get('bestsumorder_phylum', [])
        
        # Define the condition for filtering based on the specified phylum
        if phylum_filter:
            filter_condition = data_df['bestsumorder_phylum'].isin(phylum_filter)
        else:
            filter_condition = pd.Series([True] * len(data_df))  # If no filter, keep all

        # Apply the filter condition
        filtered_df = data_df[filter_condition]

        # Get the identifiers from the filtered DataFrame
        filtered_identifiers = filtered_df['identifiers'].tolist()

        # Write filtered identifiers to a text file
        with open(output_file, 'w') as file:
            file.write('\n'.join(filtered_identifiers))

        return f"Filtered identifiers saved to {output_file}"

    except Exception as e:
        return f"An error occurred: {e}"

# Arguments
data = sys.argv[1]
output_file = sys.argv[2]
taxa_filter = json.loads(sys.argv[3])  # Deserialize the taxa filter dictionary

# Call the function
status_msg = filter_identifiers(data, output_file, taxa_filter)
print(status_msg)
