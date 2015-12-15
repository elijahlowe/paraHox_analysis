from sys import argv
annot_file=open('../data/sp_gene_names.tsv')
SPU_ID={}
de={}

print "SPU ID\tWHL ID\tGene symbol\tGene description\tsuperfamily\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj"
de_file=open(argv[1])
for line in de_file:
    if line.split(',')[0]=='':
        continue
    else:
        de[line.split(',')[0]]=line.rstrip('\n').split(',')[1:]

for line in annot_file:
    key=line.split('\t')[1].rstrip(' ')
    key2=line.split('\t')[0].rstrip(' ')
    if key in de:
        print line.rstrip('\n'),'\t','\t'.join(str(p) for p in de[key])
    elif key2 in de:
        print line.rstrip('\n'),'\t','\t'.join(str(p) for p in de[key2])

