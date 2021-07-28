import cProfile, pstats, io

def clingen(infile):
	f=open(infile)
	dic={}
	for i in f:
		line=i.split("\t")
		if line[0] not in dic:
			dic[line[0]]=[[line[1]], [line[3]], [line[4]], [line[5].strip()]]
		else:
			dic[line[0]][0].append(line[1])
			dic[line[0]][1].append(line[3])
			dic[line[0]][2].append(line[4])
			dic[line[0]][3].append(line[5])
	f.close()
	return dic

def get_dd2p(infile):
	dd2p={}
	with open(infile, 'r', encoding="utf-8") as f:
		i=f.readline()
		line=i.split(",")
		if line[0] not in dd2p:
			dd2p[line[0]]=[[],[],[],[]]
		dd2p[line[0]][0].append(line[4])
		dd2p[line[0]][1].append(line[3].strip('"'))
		dd2p[line[0]][2].append(line[2].strip('"'))	
		if len(line[9])>5:
			ch=line[9].split(";")
			el=ch[:2]
			if len(el)>0:
				dd2p[line[0]][3].append('=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/'+el[0]+'","SUBSTITUTE")')
			else:
				dd2p[line[0]][3].append('na')
		else:
			dd2p[line[0]][3].append('na')
	f.close()
	return dd2p
	
def read_panel_list(infile):
	f=open(infile)
	panel=set()
	for i in f:
		line=i.strip()
		panel.add(line)
	f.close()
	return panel

