import sys
import pandas as pd
from functools import reduce
import yaml


config_file = "/home/aysevil/MolGen/faststorage/aysevil/greenland_shark/FLTransAnnot/configs/config_args.yaml"

def filter_identifiers(data, output_file, config_file):
    """Filter identifiers based on filter parameters and save to a text file.

    Args:
        data (dict): Dictionary containing data.
        output_file (str): Path to the output text file.
        filter_parameters (str): Stringified dictionary representing filter parameters.

    Returns:
        str: Message indicating the status of the operation.
    """

    #config_file="../config/config.yaml"

    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)

    print(config["blobtools_filter"])


    try:
        data_df = pd.read_csv(data, sep='\t')

        # Parse filter parameters
        #filter_params = eval(filter_parameters)
        filter_params = config["blobtools_filter"]

        # Define the order of keys
        keys_order = ['bestsumorder_superkingdom', 'bestsumorder_kingdom', 'bestsumorder_phylum', 
                      'bestsumorder_family', 'bestsumorder_class', 'bestsumorder_species']

        # Apply filters based on filter parameters
        filter_conditions = []
        for key in keys_order:
            if key in filter_params:
                filter_conditions.append(data_df[key].isin(filter_params[key]))

        # Combine filter conditions using AND operation
        if filter_conditions:
            filtered_df = data_df[reduce(lambda x, y: x & y, filter_conditions)]
        else:
            filtered_df = data_df

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

# Call the function
status_msg = filter_identifiers(data, output_file, config_file=config_file)
print(status_msg)
