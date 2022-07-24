#不能被3整除的序列(RNA gene)被过滤掉
import os
import glob
import argparse
from Bio import SeqIO


def get_orf_gc_stats(orf_seq):
    orf_codon_list = [orf_seq[i:i + 3] for i in range(0, len(orf_seq), 3)]
    total_codon_num = len(orf_codon_list)
    codon_base_1st_GC = 0
    codon_base_2nd_GC = 0
    codon_base_3rd_GC = 0
    for each_codon in orf_codon_list:
        if each_codon[0].upper() in ['G', 'C']:
            codon_base_1st_GC += 1
        if each_codon[1].upper() in ['G', 'C']:
            codon_base_2nd_GC += 1
        if each_codon[2].upper() in ['G', 'C']:
            codon_base_3rd_GC += 1

    return total_codon_num, codon_base_1st_GC, codon_base_2nd_GC, codon_base_3rd_GC


parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, help='folder name')
parser.add_argument('-x', required=True, help='file extension')
parser.add_argument('-o', required=True, help='output txt file')

args = vars(parser.parse_args())
ffn_folder  = args['i']
ffn_ext     = args['x']
output_file = args['o']


if ffn_folder[-1] == '/':
    ffn_folder = ffn_folder[:-1]

ffn_file_re = '%s/*.%s' % (ffn_folder, ffn_ext)

output_file_handle = open(output_file, 'w')

output_file_handle.write('Genome\tGC1\tGC2\tGC3\n')
for each_ffn_file in sorted(glob.glob(ffn_file_re)):

    codon_num_total = 0
    gc_1st_total = 0
    gc_2nd_total = 0
    gc_3rd_total = 0
    for each_gene in SeqIO.parse(each_ffn_file, 'fasta'):
        gene_seq = str(each_gene.seq)
        if len(gene_seq) % 3 == 0:
            current_codon_num, gc_1st, gc_2nd, gc_3rd = get_orf_gc_stats(gene_seq)
            codon_num_total += current_codon_num
            gc_1st_total += gc_1st
            gc_2nd_total += gc_2nd
            gc_3rd_total += gc_3rd
        else:
            print('skipped %s\t%sbp\t%s' % (each_ffn_file, len(gene_seq), each_gene.description))

    gc_1st_pct = float("{0:.2f}".format(gc_1st_total*100/codon_num_total))
    gc_2nd_pct = float("{0:.2f}".format(gc_2nd_total*100/codon_num_total))
    gc_3rd_pct = float("{0:.2f}".format(gc_3rd_total*100/codon_num_total))

    output_file_handle.write('%s\t%s\t%s\t%s\n' % (os.path.basename(each_ffn_file), gc_1st_pct, gc_2nd_pct, gc_3rd_pct))
output_file_handle.close()
