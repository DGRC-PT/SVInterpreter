#prepare flybase file

from sys import argv

def read_infile(infile):
	f=open(infile)
	dic={}
	for i in f:
		if i.startswith("#")==False:
			line=i.split("\t")
			if line[0] not in dic:
				dic[line[0]]=[line[3]+" "+line[5]]
			else:
				if line[3]+" "+line[5] not in dic[line[0]]:
					dic[line[0]].append(line[3]+" "+line[5])
	f.close()
	return dic

def writer(dic):
	for key,value in dic.items():
		print(key+"\t"+";".join(value))

writer(read_infile(argv[1]))
