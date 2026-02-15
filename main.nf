// include subworkflows in .nf file
include {PGS_setup} from './workflows/prscsx.nf'
include {PGS} from './workflows/prscsx.nf'

// run PRS workflow
workflow {
    prs_input = PGS_setup()
    indiv_polygenic_scores = PGS(prs_input)
}