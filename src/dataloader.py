import os

import pandas as pd


def load_natural_medicine_data(datasets=["TCM", "IMPPAT"]):
    """
    Load selected traditional medicine compound data from CSV files.
    
    Parameters:
        datasets (list): List of dataset names to include (e.g., ["TCM", "IMPPAT"])
    
    Returns:
        pd.DataFrame: Combined dataframe with a 'source' column and merged on the ID column.
    """
    data_dir = "../data/"
    available_sources = {
        "TCM": "tcm_smiles.csv",
        "IMPPAT": "imppat_smiles.csv",
    }
    
    dataframes = []
    
    for source in datasets:
        if source.upper() in available_sources:
            file_path = os.path.join(data_dir, available_sources[source.upper()])
            try:
                df = pd.read_csv(file_path)
                df["source"] = source.upper()
                df.columns.values[0] = "ID"  # Rename first column to 'ID'
                dataframes.append(df)
            except FileNotFoundError:
                print(f"File not found for source: {source}")
    
    # Merge on 'ID' if multiple dataframes, otherwise just return the single one
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df
    else:
        return pd.DataFrame()  # Return empty if nothing loaded
