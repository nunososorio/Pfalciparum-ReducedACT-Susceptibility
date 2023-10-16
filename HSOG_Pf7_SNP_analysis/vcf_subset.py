import pandas as pd
import numpy as np
import allel

# Load VCF file
callset = allel.read_vcf('hsog3.vcf')

# Convert VCF to DataFrame
hsog3_calls = allel.vcf_to_dataframe('hsog3.vcf', fields='*', alt_number=1)

# Ensure 'CHROM' is string type
hsog3_calls['CHROM'] = hsog3_calls['CHROM'].astype(str)

# Load the DR_genes.csv file
dr_genes = pd.read_csv('DR_genes.csv')

# Ensure 'Chromosome_number' is string type
dr_genes['Chromosome_number'] = dr_genes['Chromosome_number'].astype(str)

# Initialize an empty DataFrame to store the results
df_results = pd.DataFrame()

# Iterate over each row in the DR_genes DataFrame
for index, row in dr_genes.iterrows():
    # Extract the required information
    chromosome_number = row['Chromosome_number']
    start_position = row['start_position']
    end_position = row['end_position']
    
    # Filter the calls DataFrame for variants in the specified chromosome and position range
    df_filtered = hsog3_calls[(hsog3_calls['CHROM'] == chromosome_number) & 
                           (hsog3_calls['POS'] >= start_position) & 
                           (hsog3_calls['POS'] <= end_position)]
    
    # Append the filtered DataFrame to the results DataFrame
    df_results = pd.concat([df_results, df_filtered])

# Save the results to a new CSV file
df_results.to_csv('variants.csv', index=False)


