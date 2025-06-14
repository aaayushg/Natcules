import os

import pandas as pd


# Define the data_loader module
def load_natural_medicine_data():
    """
    Load traditional medicine compound data from CSV files in the main directory.
    Returns a combined DataFrame with a 'source' column to indicate the origin.
    """
    # Define file paths
    tcm_path = os.path.join("..", "tcm_smiles.csv")
    imppat_path = os.path.join("..", "imppat_smiles.csv")
    
    # Load datasets
    try:
        tcm_df = pd.read_csv(tcm_path)
        tcm_df["source"] = "TCM"
    except FileNotFoundError:
        tcm_df = pd.DataFrame()
    
    try:
        imppat_df = pd.read_csv(imppat_path)
        imppat_df["source"] = "IMPPAT"
    except FileNotFoundError:
        imppat_df = pd.DataFrame()
    
    # Combine the datasets
    combined_df = pd.concat([tcm_df, imppat_df], ignore_index=True)
    
    return combined_df

# Load data to show output
# combined_data = load_natural_medicine_data()
# import ace_tools as tools; tools.display_dataframe_to_user(name="Combined Traditional Medicine Data", dataframe=combined_data)

