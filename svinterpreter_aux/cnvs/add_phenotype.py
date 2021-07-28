from sys import argv
#hg38#GAIN
def get_categories(infile):
	f=open(infile)
	dic={}
	for i in f:
		aa=""
		line=i.split("\t")
		if float(line[10])<0.02:
			aa+="Signature; "
		if float(line[11])<0.02:
			aa+="Neuro; "
		if float(line[12])<0.02:
			aa+="Craniofacial; "
		if float(line[13])<0.02:
			aa+="Cardio; "
		if float(line[14])<0.02:
			aa+="Autism; "
		if float(line[15])<0.02:
			aa+="Epilepsy; "
		aa.strip()
		dic[line[9]]=[line[3], int(line[1]), int(line[2]), aa[:-2]]
	f.close()
	return dic

def get_Overlap(a,b):
	return max(0,min(a[1],b[1])-max(a[0],b[0]))

def complete_with_genes(infile, dic):
	f=open(infile)
	for i in f:
		gen=""
		line=i.split("\t")
		if line[4]=="copy_number_loss":
			for key, value in dic.items():
				if line[1]==value[0]:
					if get_Overlap([value[1], value[2]], [int(line[2]),int(line[3])])>=1:
						gen+=key+" ("+value[-1]+");"
			if len(gen)!=0:
				print(i.strip()+"\t"+gen[:-1])
			else:
				print(i.strip())
	f.close()

#cat table, database
complete_with_genes(argv[2], get_categories(argv[1]))
