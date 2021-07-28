#prepare wormbase file

from sys import argv


def define_DOID(infile):
	f=open(infile)
	idd=""
	dic={}
	for i in f:
		if i.startswith("id:"):
			l=i.split(" ")
			idd=l[1].strip()
		elif i.startswith("name: "):
			l=i.split(" ")
			dic[idd]=(l[1].strip())
			idd=""
	f.close()
	return dic
	
def read_infile(infile, doid):
	f=open(infile)
	dic={}
	for i in f:
		if i.startswith("!")==False and i.startswith("Taxon")==False:
			line=i.split("\t")
			if line[1]=="gene":
				if line[3] not in dic:
					dic[line[3]]=[line[8].replace("_", " ")+" "+doid[line[10]]]
				else:
					if line[8].replace("_", " ")+" "+doid[line[10]] not in dic[line[3]]:
						dic[line[3]].append(line[8].replace("_", " ")+" "+doid[line[10]])
	f.close()
	return dic

def writer(dic):
	for key,value in dic.items():
		if len(value)>0:
			print(key+"\t"+";".join(value))

writer(read_infile(argv[1], define_DOID(argv[2])))
