#! /opt/conda/bin/python

# load packages
import pandas as pd
import argparse as ap

# define arguments
def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")

    # list of all score files
    parser.add_argument('--score_files', required = True, help = 'list of all score files')
    
    return parser

args = make_arg_parser().parse_args()

# convert arguments to python variables
score_files = args.score_files

# create an empty list for dataframes to be added to
dfs = []

# reformat list of files so python can loop through it
input_file_list = list(score_files.split(" "))
print(input_file_list)

# in a loop, read in all of the PGS input files and add them to a list of dataframes
for f in input_file_list:
    print(f)
    temp_df = pd.read_table(f, sep = '\t', index_col = ['ID', 'ALLELE'])
    temp_df = temp_df[~temp_df.index.duplicated()]
    dfs.append(temp_df)

# concatenate the COLUMNS dataframes in the list (combining separate apply PGS input files)
# this makes a score file with one ID column, one allele column, and one score column for each training cohort
df = pd.concat(dfs, axis = 1)

# fill missing values with zeros (if a cohort doesn't have a SNP, we can assume the score for that SNP in that cohort is zero)
df = df.fillna(0)

# export dataframe
df.to_csv('combined_score_file.txt', sep = '\t')