import pandas as pd
from scipy.stats import norm


def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")
    
    # score file
    parser.add_argument('--scoreFile', required = True, help = 'adjusted score file')
    
    # validation population
    parser.add_argument('--valPop', required = True, help = 'validation population')
    
    # ref key
    parser.add_argument('--refKey', required = True, help = 'ref key')

    # percentile threshold
    parser.add_argument('--percentile', required = True, help = 'percentile threshold')
    
    return parser

# parse args
args = make_arg_parser().parse_args()

score_filename = args.scoreFile
val_pop = args.valPop
ref_key = args.refKey
percentile_thres = int(args.percentile)

# read in input file
score_file = pd.read_csv(score_filename, sep = '\t', low_memory = False)

# subset
score_file_sub = score_file[['IID', 'Z_norm2']]

# calculate percentiles
score_file_sub['ntile'] =  100 * norm.cdf(score_file_sub["Z_norm2"])

# classify risk
score_file_sub["risk_group"] = score_file_sub["ntile"].apply(lambda x: "high" if x >= percentile_thres else "low")

# export
score_file_sub.to_csv((val_pop + '.' + ref_key + '.pgs_percentile.txt'), sep = '\t', index = None)
