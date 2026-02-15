#! /opt/conda/bin/python

# load packages
import pandas as pd
import argparse as ap

# define arguments
def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")

    # score file name
    parser.add_argument('--scoreFile', required = True, help = 'combined score file')

    # reference panel file name
    parser.add_argument('--refPanelPvarFile', required = True, help = 'reference panel pvar file name')

    # target pvar file name
    parser.add_argument('--targetPvarFileList', required = True, help = 'list of pvar file names from all chromosomes')

    # validation population
    parser.add_argument('--valPop', required = True, help = 'validation population')
    
    return parser

args = make_arg_parser().parse_args()

# convert arguments to python variables
## score file
score_filename = args.scoreFile
## ref panel file
ref_panel_pvar_filename = args.refPanelPvarFile
## pvar file list
target_pvar_file_list = args.targetPvarFileList
## validation population
val_pop = args.valPop

# read in score file and ref panels
score_file = pd.read_csv(score_filename, sep = '\t', low_memory = False)
ref_panel_pvar_file = pd.read_csv(ref_panel_pvar_filename, sep = '\t', header = None, low_memory = False, comment = '#')

# read in chromosome separated pvar files and concat
## create an empty list for dataframes to be added to
dfs = []

## reformat list of files so python can loop through it
input_file_list = list(target_pvar_file_list.split(" "))
print(input_file_list)

## read in input files and append to df list
for f in input_file_list:
    print(f)
    dfs.append(pd.read_table(f, sep = '\t', header = None, low_memory = False, comment = '#'))

## concateante chromosome separated pvar files
pvar_file = pd.concat(dfs, axis = 0)

# make intersecting variants list
variant_list = score_file[['ID']]
variant_list = variant_list[variant_list['ID'].isin(target_pvar_file[1])]
variant_list = variant_list[variant_list['ID'].isin(ref_panel_pvar_file[1])]
print(len(variant_list.index))

# filter score file to intersecting variants
score_file_intersect = score_file[score_file['ID'].isin([variant_list['ID']])]

# export filtered score and ref panels and intersecting variant list (for plink file filtering)
score_file_intersect.to_csv((val_pop + '.combined_score_file.intersect.txt'), sep = '\t', index = None)
variant_list.to_csv((val_pop + '.variant_list.intersect.txt'), sep = '\t', index = None, header = None)