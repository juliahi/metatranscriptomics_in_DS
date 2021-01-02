import numpy as np
### REPLACING IDENTIFIERS TO GENE NAMES PREDICTED BY eggNOGmapper
slown = {}
for i in range(1,6):#6 - number of files 
	f = open('./eggNOGmapper/part'+str(i)+'_no_tab.faa.emapper.annotations', 'r').readlines()
	for line in f[:]:
		line = line.strip().split('\t')
		ident = (line[0].split('|')[0], line[0].split('|')[6]) 
		print (ident)
		slown [ident] = line[4] # dictionary with identifiers and gene names
		print (slown [ident] )
print (len(slown))
#print (slown)

g = open('./all4_all16_pvals.csv', 'r').readlines()
h = open('./all4_all16_pvals_genes.csv', 'w')
h.write(g[0])
counter =0
print('\n')
for ln in g[1:]:
	ln = ln.strip().split(',')
	id2 = (ln[0].split('|')[0], ln[0].split('|')[6])
	print (id2)
	# slown[id2] gene name
	if id2 in slown:
		print (slown[id2])
		if len(slown[id2])>0:
			h.write(slown[id2]+'\t'+'\t'.join(ln[1:])+'\n')
		if len(slown[id2])==0:
			h.write('.\t'+'\t'.join(ln[1:])+'\n')
	else:
		counter+=1
print ('This number of putative genes were not predicted by eggnogmapper',counter)

### NORMALIZE COUNTS
l = open('./metagenemark_counts_filtered_genes_all.txt', 'r').readlines()
mice = l[0].strip().split('\t')[1:]
k = open('./metagenemark_counts_filtered_genes_all_normalised.txt', 'w')
k.write(l[0])
all_quartiles = []
all_exprs =[]
genes = [ ln.strip().split('\t')[0] for ln in l[1:]]

for i in range(1, len(mice)+1):
	exprs = []
	for line in l[1:]:
		line=line.strip().split('\t')
		exprs.append(float(line[i]))
		#print (line[0])

	quartile = np.percentile(exprs, 75)
	all_quartiles.append(quartile)
	all_exprs.append(exprs)
quart = np.median(all_quartiles)

print(len(exprs))
print(len(genes))


all_exprs_norm = []
for j in range(len(mice)):
	quart_lokal = all_quartiles[j]
	exprs_norm =[]
	mouse = all_exprs[j]
	for ex in mouse:
		ex_norm = ex*quart/quart_lokal*1.0
		exprs_norm.append(ex_norm)
	all_exprs_norm.append(exprs_norm)

for nr in range(len(genes)):
	gene = genes[nr]
	k.write(gene+'\t')
	what_to_write = [str(all_exprs_norm[o][nr]) for o in range(len(mice))]
	k.write('\t'.join(what_to_write)+'\n')



### FILTER GENES WITH KNOWN SYMBOLS AND AVERAGE VALUES FOR GENES WITH MULTIPLE OCCURRENCES
l = open('./metagenemark_counts_filtered_genes_all_normalised.txt', 'r').readlines()
mice = l[0].strip().split('\t')[1:]
genes = set([ln.strip().split('\t')[0] for ln in l[1:] if ln.strip().split('\t')[0]!='.']) # wybieram nazwy genow, oprocz tych ktore nie mialy nazw (mialy kropki)

slown2 = {} # dictionary for keeping the calculated averages
for mouse in mice:
	for gene in genes:
		slown2[(mouse, gene)] = []
print(len(mice))	

for ln in l[1:]:
	ln = ln.strip().split('\t')
	gene = ln[0]
	if gene!='.':
		#print (ln)
		for i in range(1,len(mice)+1):
			#print (i)
			mouse = l[0].strip().split('\t')[i]	
			#print (mouse, gene)
			slown2[(mouse, gene)].append(float(ln[i]))
slown_average = {}
for el in slown2:
	slown_average[el] = np.mean(slown2[el])

m = open('./metagenemark_counts_filtered_genes_all_normalised_averaged.txt', 'w')

m.write('gene\t'+'\t'.join(mice)+'\n')
for gene in genes:
	what_to_write = [slown_average[(mouse, gene)] for mouse in mice]
	what_to_write_str = [str(el1) for el1 in what_to_write]
	m.write(gene+'\t'+'\t'.join(what_to_write_str)+'\n')'''

			

