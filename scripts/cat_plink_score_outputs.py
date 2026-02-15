# load packages
import pandas as pd
import os
import argparse as ap
import sys

# create arguments
def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")
    
    parser.add_argument('-v', '--valPop', required = True, help = 'Validation population')
    
    parser.add_argument('--score_file_list', required = True, help = 'List of chromosome separated score files for each validation population')

    parser.add_argument('--variant_list', required = True, help = 'List of chromosome separated variant lists used in score computation for each validation population')

    return parser

args = make_arg_parser().parse_args()

# define arguments as variables
validation_population = args.valPop
score_file_list = args.score_file_list
variant_list = args.variant_list

# reformat list of files that will work in a for loop- score files and variant lists
# score files
score_file_string = ','.join(score_file_list)
score_file_string = score_file_string.replace('[', '')
score_file_string = score_file_string.replace(']', '')
score_file_string = score_file_string.replace(',', '')
score_file_list_reformatted = list(score_file_string.split(" "))
print(score_file_list_reformatted)

# variant lists
variant_list_string = ','.join(variant_list)
variant_list_string = variant_list_string.replace('[', '')
variant_list_string = variant_list_string.replace(']', '')
variant_list_string = variant_list_string.replace(',', '')
variant_list_reformatted = list(variant_list_string.split(" "))
print(variant_list_reformatted)

# concatenate variant lists
print(f"concatenating lists of variants using in plink --score for {validation_population}")

# create empty list of dataframes to append variant lists to
variant_dfs = []

# read in variant lists in a for loop and append to a list of dataframes
for f in variant_list_reformatted:
    variant_dfs.append(pd.read_table(f, sep = '\t', header = None))

# concatenate dataframes in the list (and therefore variant lists for all the chromosomes)
variant_list_cat = pd.concat(variant_dfs)

# sort concatenated dataframe
variant_list_cat.sort_values(by = 0, inplace = True)

# add a column name
variant_list_cat.rename(columns = {0 : 'Variant_ID'}, inplace = True)

# export dataframes
variant_list_cat.to_csv(f'{validation_population}.all_computed_PGS.variant_list.txt', sep = '\t', index = False)


# concatenate plink score outptpus
print(f"concatenating plink score output files for {validation_population}")

# create empty score dataframe
score_cat_df = None

# read in score files in a for loop and concatenate
for f in score_file_list_reformatted:
    score_df = pd.read_table(f)
    if '#FID' in score_df.columns and '#IID' in score_df.columns:
        score_df.set_index(['#FID', '#IID'], drop = True, inplace = True)
    elif '#FID' in score_df.columns and 'IID' in score_df.columns:
        score_df.set_index(['#FID', 'IID'], drop = True, inplace = True)
    elif '#IID' in score_df.columns:
        score_df.set_index(['#IID'], drop = True, inplace = True)
    else:
        raise ValueError('ID columns do not meet conditions')
    # take the sum of the columns in all the dataframes- thereby adding up the chromosome separated scores
    all_cols = score_df.columns
    if score_cat_df is None:
        score_cat_df = score_df[all_cols]
    else:
        score_cat_df = score_cat_df.add(score_df[all_cols], fill_value = 0)

# create an average score for each person by dividing the score sum / allele count
sum_cols = [c for c in score_cat_df.columns if 'SCORE_SUM' in c]
avg_cols = []
for c in sum_cols:
    c = c.replace('SCORE_SUM', 'SCORE_AVG')
    print(c)
    avg_cols.append(c)
score_cat_df[avg_cols] = score_cat_df[sum_cols].div(score_cat_df['ALLELE_CT'], axis = 0)

# reorder columns
score_cat_df_no_order = score_cat_df[['ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM']]
score_cat_df.drop(columns = ['ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM'], inplace = True)
score_cat_df = score_cat_df.reindex(sorted(score_cat_df.columns), axis = 1)
score_cat_df_ordered = pd.concat([score_cat_df_no_order, score_cat_df], axis = 1)

# export dataframe
score_cat_df_ordered.to_csv(f'{validation_population}.all_computed_PGS.txt', sep = '\t')