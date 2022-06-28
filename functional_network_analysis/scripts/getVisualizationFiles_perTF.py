import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get files for network visualization on Cytoscape',
                        conflict_handler='resolve')

    parser.add_argument('enr_out_info', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif enrichment output file')
    
    parser.add_argument('desc_ortho', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with genes description '+\
                        'and Arabidopsis orthology information')

    parser.add_argument('mot_tf', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with motif information')

    parser.add_argument('tf_id', nargs = 1, type = str,
                    help = '',
                    metavar = 'TF ID to retrieve network')
    
    parser.add_argument('net_file', nargs = 1, type = str,
                    help = '',
                    metavar = 'Output network file')
    
    parser.add_argument('node_file', nargs = 1, type = str,
                    help = '',
                    metavar = 'Output node annotation file')
    
    parser.add_argument('-s', '--set_id', nargs = '?', type = str,
                    help = '', default = False,
                    metavar = 'If there is more than one set ' +\
                        'specify for which set to retrieve the ' +\
                        'network from')
    
    args = parser.parse_args()

    return args

args = parseArgs()

enr_file = args.enr_out_info[0]
desc_ortho_file = args.desc_ortho[0]
mot_tf_file = args.mot_tf[0]
tf = args.tf_id[0]
net_file_out = args.net_file[0]
node_file_out = args.node_file[0]
set_id = args.set_id

enr_info = pd.read_csv(enr_file, sep = "\t")

ortho_desc_df = pd.read_csv(desc_ortho_file, sep = "\t")

tfs_with_mot = set(pd.read_csv(mot_tf_file, sep = "\t")['gene_id'].unique())

assert tf in set(enr_info.SolycID.unique()), "TF not enriched, network cannot be generated"

if set_id == False:
    tgs = set(enr_info[enr_info.SolycID == tf]['Target genes'].str.split(",", expand = True).melt().dropna()['value'])
    
    sets = set(enr_info['Set ID'].unique())
    
    assert len(sets) == 1, "Please specify for which set you want to generate the network"
    
else:
    tgs = set(enr_info[(enr_info.SolycID == tf) & (enr_info['Set ID'] == set_id)]['Target genes'].str.split(",", expand = True).melt().dropna()['value'])
    
with open(net_file_out, 'w') as fout:
    for tg in tgs:
        fout.write("\t".join([tf, tg]))
        fout.write("\n")
        
all_genes = tgs.union(set([tf]))
nodes_df = ortho_desc_df[ortho_desc_df.SolycID.isin(all_genes)]
nodes_df.loc[:, 'is_regulator'] = (nodes_df.SolycID == tf)
nodes_df.loc[:, 'is_TF'] = (nodes_df.SolycID.isin(tfs_with_mot))
nodes_df.to_csv(node_file_out, sep = "\t", index = None)