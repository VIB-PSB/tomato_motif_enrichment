import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import argparse


def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get network from motif enrichment output',
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
    
    parser.add_argument('-f', '--filter_by_set', nargs = '?', type = str,
                    help = '', default = False,
                    metavar = 'Filter network for genes present in file')
    
    args = parser.parse_args()

    return args

args = parseArgs()

enr_file = args.enr_out_info[0]
desc_ortho_file = args.desc_ortho[0]
mot_tf_file = args.mot_tf[0]
net_file_out = args.net_file[0]
node_file_out = args.node_file[0]
set_id = args.set_id
set_filter_file = args.filter_by_set

D = {}
with open(enr_file, 'r') as fin:
    for line in fin:
        if line.startswith("Set ID"):
            continue
        rec = line.strip().split("\t")
        gene_set, TF, TGs = rec[0], rec[8], rec[-1].split(",")
        tmp_set = set()
        for TG in TGs:
            tmp_set.add(tuple([TF, TG]))
        D.setdefault(gene_set, []).append(tmp_set)

if set_filter_file != False:
    set_filter = set()
    with open(set_filter_file, 'r') as f:
        for line in f:
            rec = line.strip().split("\t")
            set_filter.add(rec[0])
        
if set_id == False:
    
    assert len(D) == 1, "Please specify for which set you want to generate the network"
    
    for i in D:
        with open(net_file_out, 'w') as fout:
            net = set.union(*D[i])
            for tf_tg_pair in net:
                if set_filter_file == False:
                    fout.write("\t".join([tf_tg_pair[0], tf_tg_pair[1]]))
                    fout.write("\n")
                else:
                    tf, tg = tf_tg_pair[0], tf_tg_pair[1]
                    if tf in set_filter and tg in set_filter:
                        fout.write("\t".join([tf_tg_pair[0], tf_tg_pair[1]]))
                        fout.write("\n")
else:
    for i in D:
        if i != set_id:
            continue
        with open(net_file_out, 'w') as fout:
            net = set.union(*D[i])
            for tf_tg_pair in net:
                if set_filter_file == False:
                    fout.write("\t".join([tf_tg_pair[0], tf_tg_pair[1]]))
                    fout.write("\n")
                else:
                    tf, tg = tf_tg_pair[0], tf_tg_pair[1]
                    if tf in set_filter and tg in set_filter:
                        fout.write("\t".join([tf_tg_pair[0], tf_tg_pair[1]]))
                        fout.write("\n")

                
ortho_desc_df = pd.read_csv(desc_ortho_file, sep = "\t")

tfs_with_mot = set(pd.read_csv(mot_tf_file, sep = "\t")['gene_id'].unique())

enr_info = pd.read_csv(enr_file, sep = "\t")

if set_id == False:
    tgs = set(enr_info['Target genes'].str.split(",", expand = True).melt().dropna()['value'])
    
    sets = set(enr_info['Set ID'].unique())
    
    assert len(sets) == 1, "Please specify for which set you want to generate the network"
    
else:
    tgs = set(enr_info[enr_info['Set ID'] == set_id]['Target genes'].str.split(",", expand = True).melt().dropna()['value'])
    
all_genes = tgs.union(set(enr_info.SolycID.unique()))
nodes_df = ortho_desc_df[ortho_desc_df.SolycID.isin(all_genes)]
nodes_df.loc[:, 'is_TF'] = (nodes_df.SolycID.isin(tfs_with_mot))
nodes_df.to_csv(node_file_out, sep = "\t", index = None)