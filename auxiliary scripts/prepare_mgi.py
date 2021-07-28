#prepare MGI file

from sys import argv

def read_MP(infile):
	f=open(infile)
	dic={}
	for i in f:
		line=i.split("\t")
		dic[line[0]]=line[1]
	f.close()
	return dic

def read_infile(infile, mps):
	f=open(infile)
	dic={}
	for i in f:
		if i.startswith("#")==False:
			line=i.split("\t")
			if line[-2] not in dic:
				dic[line[-2]]=[mps[line[4]]]
			else:
				if mps[line[4]] not in dic[line[-2]]:
					dic[line[-2]].append(mps[line[4]])
	f.close()
	return dic

def writer(dic):
	for key,value in dic.items():
		print(key+"\t"+";".join(value))

writer(read_infile(argv[1], read_MP(argv[2])))
