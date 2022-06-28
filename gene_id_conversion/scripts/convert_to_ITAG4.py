import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to convert ITAG 2.5 gene IDs to ITAG 4.0',
                        conflict_handler='resolve')

    parser.add_argument('gene_list', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with ITAG 2.5 gene IDs')
    
    parser.add_argument('conv_table', nargs = 1, type = str,
                        help = '',
                        metavar = 'Conversion table ITAG 2.5 to ITAG 4.0')
    
    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Output file')
    
    args = parser.parse_args()

    return args

args = parseArgs()

infile = args.gene_list[0]
gene_id_mapping_table = args.conv_table[0]
outfile = args.output[0]

v2_v4_mappings = {}
with open(gene_id_mapping_table, 'r') as f:
    for line in f:
        rec = line.strip().split("\t")
        v4_gene = rec[1]
        v3_gene = rec[0]
        v2_v4_mappings.setdefault(v3_gene, set()).add(v4_gene)

with open(outfile, 'w') as fout:
    with open(infile, 'r') as fin:
        for line in fin:
            rec = line.strip().split("\t")
            gene_v2 = rec[0]
            if gene_v2 in v2_v4_mappings:
                for gene_v4 in v2_v4_mappings[gene_v2]:
                    fout.write(gene_v4)
                    fout.write("\n")
