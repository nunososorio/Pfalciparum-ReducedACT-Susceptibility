# -*- coding: utf-8 -*-
import malariagen_data
import pandas as pd

pf7 = malariagen_data.Pf7()
variant_dataset = pf7.variant_calls()

with open('Pf7_sample_index_FWS_095.txt', 'r') as file:
    first_column_list = file.readlines()

first_column_list = [element.strip() for element in first_column_list]
first_column_list = [int(element) for element in first_column_list]

subset = variant_dataset.isel(samples = first_column_list)

positions = [489337, 375427, 383169, 596674, 2481108]
chroms = ['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_09_v3', 'Pf3D7_13_v3']

indexer = False
for pos, chrom in zip(positions, chroms):
    indexer |= ((variant_dataset.variant_position == pos) & (variant_dataset.variant_chrom == chrom)).compute()

subset2 = subset.isel(variants=indexer)
subset2_df = subset2.to_dataframe()

metadata = pf7.sample_metadata()

geosnps = subset2_df.merge(metadata[['Sample', 'Population']], left_on='sample_id', right_on='Sample', how='left')
geosnps = geosnps.drop_duplicates()

geosnps_gen1 = geosnps[geosnps['call_genotype'] == 1]

total_count = geosnps.groupby(['variant_position', 'Population']).size()
filtered_count = geosnps_gen1.groupby(['variant_position', 'Population']).size()

percentage_df = ((filtered_count / total_count) * 100).round(2)
percentage_df = percentage_df.reset_index().sort_values(['variant_position', 0], ascending=[True, False]).set_index(['variant_position', 'Population'])

# Save the DataFrame as a CSV file
percentage_df.to_csv('percentage.csv')