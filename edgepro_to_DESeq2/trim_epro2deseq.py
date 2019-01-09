''
parse edge-pro combined rpkm file for duplicates
usage: python trim_epro2deseq.py .csv
'''

import sys

# command line parameters
csv = sys.argv[1]

# initialize variables
gene_ids = []
count_totals = []
lines = []
firstRow = True

with open(csv, 'r') as f:
    for line in f:
        # add current line to list of lines
        lines.append(line)
        if firstRow == True:
            firstRow = False
            # remove 'gene' header such that R can import the gene ids as rows
            # remove file extensions
            lines[0] = lines[0].replace('gene', '')
            lines[0] = lines[0].replace('.out.rpkm_0', '')
        line = line.replace('\n', '')
        gene = line.split('\t')[0]
        # calculate total count for a row, sum of columns 1 to 4
        if (gene != 'gene'):
            count = int(line.split('\t')[1]) + int(line.split('\t')[2]) + int(line.split('\t')[3]) + int(line.split('\t')[4])
        # if gene is a duplicate:
        # remove old entry in lines if new count > old count
        # remove new entry in lines if old count > new count
        if (gene_ids != []) & (gene in gene_ids):
            if count > count_totals[-1]:
                lines.pop(-2)
            else:
                lines.pop()
        # keep track of gene_ids and counts
        if (gene.strip() != '') & (gene != 'gene'):
             gene_ids.append(gene)
             count_totals.append(count)

# write parsed lines to final deseqFile
o = open('Listeria_deseqFile','w')
for line in lines:
    o.write(line)
o.close()
