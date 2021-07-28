#prepare rat file

from sys import argv


def define_associatedwith(v):
	beg=[]
	for el in v.split(";"):
		if "ASSOCIATED" in el:
			beg.append(el.replace("ASSOCIATED WITH", ""))
		elif len(beg)>0:
			if "INTERACTS" not in el and "FOUND" not in el:
				beg.append(el)
			else:
				break
	return beg
	
def read_infile(infile):
	f=open(infile)
	dic={}
	for i in f:
		if i.startswith("#")==False and i.startswith("GENE_")==False:
			line=i.split("\t")
			dic[line[0]]=define_associatedwith(line[3])
	f.close()
	return dic

def writer(dic):
	for key,value in dic.items():
		if len(value)>0:
			print(key+"\t"+";".join(value))

writer(read_infile(argv[1]))
