from sys import argv
annot_file=open('../data/pm_gene_names.tsv')
SPU_ID={}
de={}

de_file=open(argv[1])

print 'PMI ID\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tGene symbol'
for line in de_file:
    de[line.split(',')[0]]=line.rstrip('\n').split(',')[1:]

for line in annot_file:
    key=line.split('\t')[0].rstrip(' ')
    key2=line.split('\t')[1].rstrip(' ')
    if key in de:
        print line.rstrip('\n').split('\t')[0]+'\t'+'\t'.join(str(p) for p in de[key])+'\t'+ line.rstrip('\n').split('\t')[1]
    elif key2 in de:
        print line.rstrip('\n')+'\t'+'\t'.join(str(p) for p in de[key2])

