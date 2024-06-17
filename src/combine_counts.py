
import pandas as pd
import glob

# List of files to combine
files = sorted(glob.glob("barcode*_filtered_rel-abundance.tsv"))

# Initialize an empty DataFrame to store the combined data
combined_df = pd.DataFrame()

# Define the metadata columns
metadata_cols = ["tax_id", "species", "genus", "family", "order", "class", "phylum", "clade", "superkingdom", "subspecies",$

# Iterate through each file
for file in files:
    # Extract the barcode from the filename
    barcode = file.split("_")[0]

    # Read the current TSV file
    df = pd.read_csv(file, sep="\t")

    # Rename the "estimated counts" column to the barcode
    df = df.rename(columns={"estimated counts": barcode})

    # If combined_df is empty, initialize it with the first DataFrame
    if combined_df.empty:
        combined_df = df
    else:
        # Merge the current DataFrame with the combined DataFrame on the "tax_id" column
        combined_df = pd.merge(combined_df, df, on=["tax_id"], how="outer", suffixes=('', '_dup'))

        # Resolve duplicated columns by keeping only one of the duplicates
        combined_df = combined_df.loc[:,~combined_df.columns.duplicated()]

# Ensure all metadata columns are included in the combined DataFrame
for col in metadata_cols:
    if col not in combined_df.columns:
        combined_df[col] = None

# Reorder columns: metadata columns first, followed by barcode columns
barcode_cols = [col for col in combined_df.columns if col not in metadata_cols]
combined_df = combined_df[metadata_cols + barcode_cols]

# Save the combined DataFrame to a CSV file
combined_df.to_csv("combined_estimated_counts.csv", index=False)

