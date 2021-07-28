#prepare zfin file

from sys import argv

def read_infile(infile):
	f=open(infile)
	dic={}
	f.readline()
	f.readline()
	for i in f:
		line=i.split("\t")
		if line[0] not in dic:
			if len(line[4])>2:
				dic[line[0]]=[line[4]]
			else:
				dic[line[0]]=[line[6]]
		else:
			if line[4] not in dic[line[0]] and line[6] not in dic[line[0]]:
				if len(line[4])>2:
					dic[line[0]].append(line[4])
				else:
					dic[line[0]].append(line[6])
	f.close()
	return dic

def writer(dic):
	for key,value in dic.items():
		print(key+"\t"+";".join(value))

writer(read_infile(argv[1]))
