#! /opt/conda/bin/python

# load packages
import pandas as pd
import argparse as ap
import sys
import os

# define arguments
def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")

    # pgs id
    parser.add_argument('--pgsId', required = True, help = 'PGS weight ID')

    # score file name
    parser.add_argument('--scoreFile', required = True,  help = 'score file')

    # chromosome column name in score file
    parser.add_argument('--scoreChrCol', required = True, help = 'position/base pair column name in score files')

    # position column name in score file
    parser.add_argument('--scorePosCol', required = True, help = 'position/base pair column name in score files')

    # A1 columns name in score file
    parser.add_argument('--scoreA1Col', required = True, help = 'A1 column name in score files')

    # A2 column name in score file
    parser.add_argument('--scoreA2Col', required = True, help = 'A2 column name in score files')
    
    # PGS column name in score file
    parser.add_argument('--scorePGSCol', required = True, help = 'PGS column name in score files')
    
    return parser


args = make_arg_parser().parse_args()

# convert arguments to python variables
## pgs id
pgs_id = args.pgsId
## chromosome
chr_colname = args.scoreChrCol
## position
pos_colname = args.scorePosCol
## A1
a1_colname = args.scoreA1Col
## A2
a2_colname = args.scoreA2Col
## PGS
pgs_colname = args.scorePGSCol
## score file name
score_file_name = args.scoreFile

# add score variant ID to list
var_id_type_list.append(args.scoreIdFormat)

# read in score file
score_file = pd.read_table(score_file_name, sep = None, engine = 'python', dtype = {pos_colname: int, variant_id_colname: str, a1_colname: str, a2_colname: str})
    
# get column order variable
og_col_order = score_file.columns

# reformat variant ID
score_file['new_variantID'] = 'chr' + score_file['chr_colname'].astype(str) + ':' + score_file['pos_colname'].astype(str) + ':' score_file['a2_colname'].astype(str) + ':' score_file['a1_colname'].astype(str)

# subset
score_file_sub = score_file[['new_variantID', a1_colname, pgs_colname]]

# rename columns
score_file_sub = new_score_file_sub.rename(columns = {'new_variantID' : 'ID',
                                                        a1_colname : 'ALLELE',
                                                        pgs_colname : f'{pgs_id}.SCORE'})

# export file
score_file_sub.to_csv(f'${pgs_id}.reformatted.txt', sep = '\t', index = None)
    
    