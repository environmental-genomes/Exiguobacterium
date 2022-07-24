from Bio import SeqIO
import sys

if len(sys.argv) != 2:
    print("Usage: need a assembly file")
    exit()

assembly = sys.argv[1]

i = 0
for record in SeqIO.parse(assembly, format="fasta"):
    i = i + 1
    record.id = assembly[0:14] + 'c' + str(i)
    print('>' + record.id)
    print(record.seq)
