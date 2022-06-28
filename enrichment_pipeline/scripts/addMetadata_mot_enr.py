import pandas as pd
import argparse


def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to add metadata to TF2Network output',
                        conflict_handler='resolve')

    parser.add_argument('enr_out', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif enrichment output file')
    
    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Name of the output file')
    
    parser.add_argument('metadata', nargs = 1, type = str,
                    help = '',
                    metavar = 'file with tomato motif metadata ')

    args = parser.parse_args()

    return args

args = parseArgs()

enr_file = args.enr_out[0]
outfile = args.output[0]
metadata_file = args.metadata[0]

enr_df = pd.read_csv(enr_file, sep = "\t", comment = "#", header = None)
enr_df.columns = ['set_id', 'ftr_id', 'p_val', 'q_val', 'enr_fold', 'set_size', 'ftr_size', 'n_hits', 'hits']

metadata_df = pd.read_csv(metadata_file, sep = "\t")

enr_df_info = enr_df.merge(metadata_df, how = 'left', left_on = 'ftr_id', right_on = 'motif_id')
enr_df_info["enr_rank"] = enr_df_info.groupby('set_id')["q_val"].rank("dense")
enr_df_info = enr_df_info.drop('motif_id', axis = 1)

enr_df_info.columns = ['Set ID', 'Motif ID', 'p-value', 'q-value', 'Enrichment fold', 'Number of genes of interest', 'Number of genes targeted by motif', 'Number of hits', 'Target genes', 'SolycID', 'TF description', 'TF/motif family', 'Arabidopsis ortholog gene ID', 'Arabidopsis ortholog gene symbol', 'Enrichment rank']

enr_df_info = enr_df_info[['Set ID', 'Motif ID', 'p-value', 'q-value', 'Enrichment fold', 'Number of genes of interest', 'Number of genes targeted by motif', 'Number of hits', 'SolycID', 'TF description', 'TF/motif family', 'Arabidopsis ortholog gene ID', 'Arabidopsis ortholog gene symbol', 'Enrichment rank', 'Target genes']]

enr_df_info.to_csv(outfile, sep = "\t", index = None)