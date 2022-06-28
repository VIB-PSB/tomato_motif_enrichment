import argparse


def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get network from TF2Network output',
                        conflict_handler='resolve')

    parser.add_argument('enr_out_info', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif enrichment output file')
    
    parser.add_argument('tf_id', nargs = 1, type = str,
                        help = '',
                        metavar = 'TF ID to retrieve TGs')
    
    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'output file')
    
    args = parser.parse_args()

    return args

args = parseArgs()

enr_info = args.enr_out_info[0]
desired_tf = args.tf_id[0]
out_file = args.output[0]

tgs_set = set()

with open(enr_info, 'r') as fin:
    for line in fin:
        rec = line.strip().split("\t")
        tf_id = rec[8]
        tgs = set(rec[14].split(","))
        if tf_id == desired_tf:
            tgs_set = tgs_set.union(tgs)
            
with open(out_file, 'w') as fout:
    for gene in tgs_set:
        fout.write(gene)
        fout.write("\n")