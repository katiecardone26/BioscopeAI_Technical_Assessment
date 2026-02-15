#! /opt/conda/bin/python

# load packages
import pandas as pd
import argparse as ap

# define arguments
def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")

    # sample list
    parser.add_argument('--sample_list', required = True, help = 'sample list file')

    # psam file
    parser.add_argument('--psam_file', required = True, help = 'one psam file')

    # sample file delimiter
    parser.add_argument('--sample_delim', required = True, help = 'delimiter in sample file')

    # sample id col
    parser.add_argument('--sample_id_col', required = True, help = 'id col name in sample file')

    # sample id col
    parser.add_argument('--valPop', required = True, help = 'validation population')

    return parser

args = make_arg_parser().parse_args()

# convert arguments to python variables
## sample list filename
sample_list_filename = args.sample_list
## psam filename
psam_filename = args.psam_file
## sample file delimiter
sample_delim = args.delim
## sample id col
sample_id_col = args.sample_id_col
## val pop
val_pop = args.valPop

# import psam file
psam_file = pd.read_csv(psam_filename, sep = None, engine = 'python', header = None, comment = '#', usecols = [0, 1], names = ['FID', 'IID'], dtype = str)

# reformat psam file to just IID
if psam['IID'].isnull().all():
    psam = psam.drop(columns = ['IID'])
    psam['IID'] = psam['FID']
    psam = psam.drop(columns = ['FID'])

# check if sample file has a header
if sample_id_col.isdigit():
    # sample file does not have a header
    ## read in file
    sample = pd.read_csv(sample_list_filename, sep = sample_delim, engine = 'python', header = None, dtype = str)
    ## filter psam file to IDs in sample file
    psam_sample = psam[psam['IID'].isin(sample[sample_id_col])]
else:
    # sample file has a header
    ## read in file
    sample = pd.read_csv(sample_list_filename, sep = sample_delim, engine = 'python', dtype = str)
    ## filter fam/psam file to IDs in sample file
    psam_sample = psam[psam['IID'].isin(sample[sample_id_col])]

# export sample list
psam_sample.to_csv((val_pop + '.fam_format_sample_list.txt'), sep = '\t', index = None, header = None)