from pybedtools import *
import os
import argparse


def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to process process liftoff output ' + \
                                            'and intersect with query genome GFF and '+ \
                                            'generate the final gene ID mapping file.',
                        conflict_handler='resolve')

    parser.add_argument('liftoff_output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Liftoff GFF output file')
    
    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Name of the output file')
    
    parser.add_argument('target_gff_filtered', nargs = 1, type = str,
                    help = '',
                    metavar = 'GFF of the target genome assembly '+ \
                              'filtered for gene features')

    parser.add_argument('-cov', '--coverage_cutoff', nargs = '?', type = float,
                        default = 0.75,
                        metavar = 'Coverage cutoff for liftoff alignments')
    
    parser.add_argument('-ov', '--bedtools_overlap', nargs = '?', type = float,
                        default = 0.90,
                        metavar = 'bedtools overlap cutoff for intersection')

    args = parser.parse_args()

    return args

args = parseArgs()

liftoff_gff = args.liftoff_output[0]
outfile = args.output[0]
target_gff_genes = args.target_gff_filtered[0]
cutoff = args.coverage_cutoff
cov = args.bedtools_overlap

tmp_file = 'tmp_liftoffOutput_filtered.txt'

with open(tmp_file, 'w') as fout:
    with open(liftoff_gff, 'r') as fin:
        for line in fin:
            rec = line.strip().split("\t")
            if rec[2] == 'gene':
                last_col = rec[-1].split(";")
                coverage = float(last_col[5].split("=")[1])
                if coverage >= cutoff:
                    fout.write(line)
                    
itag4_gff = BedTool(target_gff_genes)

bed_file = BedTool(tmp_file)
D = {}

int1 = bed_file.intersect(itag4_gff, f = cov, wb = True, wa = True)
int2 = bed_file.intersect(itag4_gff, F = cov, wb = True, wa = True)
with open(int1.fn, 'r') as fint1:
    for line in fint1:
        rec = line.strip().split("\t")
        gene_itag2 = dict(item.split("=") for item in rec[8].split(";"))['Name']
        gene_itag4 = dict(item.split("=") for item in rec[-1].split(";"))['Name']
        if gene_itag2 not in D:
            D[gene_itag2] = set([gene_itag4])
        elif gene_itag2 in D:
            if gene_itag4 not in D[gene_itag2]:
                D[gene_itag2].add(gene_itag4)
            else:
                pass
            
with open(int2.fn, 'r') as fint2:
    for line in fint2:
        rec = line.strip().split("\t")
        gene_itag2 = dict(item.split("=") for item in rec[8].split(";"))['Name']
        gene_itag4 = dict(item.split("=") for item in rec[-1].split(";"))['Name']
        if gene_itag2 not in D:
            D[gene_itag2] = set([gene_itag4])
        elif gene_itag2 in D:
            if gene_itag4 not in D[gene_itag2]:
                D[gene_itag2].add(gene_itag4)
            else:
                pass

D_4 = {}

with open(int1.fn, 'r') as fint1:
    for line in fint1:
        rec = line.strip().split("\t")
        gene_itag2 = dict(item.split("=") for item in rec[8].split(";"))['Name']
        gene_itag4 = dict(item.split("=") for item in rec[-1].split(";"))['Name']
        if gene_itag4 not in D_4:
            D_4[gene_itag4] = set([gene_itag2])
        elif gene_itag4 in D_4:
            if gene_itag2 not in D_4[gene_itag4]:
                D_4[gene_itag4].add(gene_itag2)
            else:
                pass
            
with open(int2.fn, 'r') as fint2:
    for line in fint2:
        rec = line.strip().split("\t")
        gene_itag2 = dict(item.split("=") for item in rec[8].split(";"))['Name']
        gene_itag4 = dict(item.split("=") for item in rec[-1].split(";"))['Name']
        if gene_itag4 not in D_4:
            D_4[gene_itag4] = set([gene_itag2])
        elif gene_itag4 in D_4:
            if gene_itag2 not in D_4[gene_itag4]:
                D_4[gene_itag4].add(gene_itag2)
            else:
                pass

new_D = dict()
for gene_map in D:
    if len(D[gene_map]) == 1:
        new_D[gene_map, list(D[gene_map])[0]] = '1-to-1'
    else:
        for gene4 in D[gene_map]:
            new_D[gene_map, gene4] = '1-to-many'


new_D_4 = dict()
for gene_map in D_4:
    if len(D_4[gene_map]) == 1:
        new_D_4[list(D_4[gene_map])[0], gene_map] = '1-to-1'
    else:
        for gene4 in D_4[gene_map]:
            new_D_4[gene4, gene_map] = 'many-to-1'

with open(outfile, 'w') as fout:
    for mapping in new_D:
        
        norm = new_D[mapping]
        inv = new_D_4[mapping]
        
        if 'many' in norm and 'many' in inv:
            rel = 'many-to-many'
        elif 'many' in norm and 'many' not in inv:
            rel = norm
        elif 'many' in inv and 'many' not in norm:
            rel = inv
        else:
            rel = '1-to-1'
        
        fout.write("\t".join([mapping[0], mapping[1], rel]))
        fout.write("\n")
        
os.remove("tmp_liftoffOutput_filtered.txt")