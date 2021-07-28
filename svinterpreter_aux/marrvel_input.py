import cProfile, pstats, io
import requests
from sys import argv
import re

def get_description(gene_id):
	"""recives the gene_id and retrives the description
	of the function. If the text has less than 281 characters
	is fully written. Otherwise the number of sentences is shortened
	from the end to the begining, until it attains less than 281 characteres.
	returns the function text."""
	url = "http://v1.marrvel.org/data/OMIM/"
	req = requests.get(url, params = {"geneSymbol": gene_id})
	data=req.json()
	gg=""
	if data!=None and "description" in data:
		temp=re.sub("[\(\[].*?[\)\]]", "",data["description"])
		if len(temp)<=280:
			gg=temp
		else:
			aa=temp.split(".")
			b=len(aa)-2
			res=""
			while b>0:
				if len(".".join(aa[:b]))<=280:
					res=".".join(aa[:b])
					break
				b-=1
			if res=="":
				res=aa[0]
			gg=res+"."
	return gg

def get_omim(gene_id):
	inh={"Autosomal dominant":"AD", "Autosomal recessive":"AR", "Isolated cases":"IC","Somatic mutation":"SMu", "Multifactorial":"Mu", "Mitochondrial": "Mi", "Somatic mosaicism":"SMo", "X-linked":"XL","X-linked dominant":"XLD", "X-linked recessive":"XLR", "Digenic recessive":"DR", "Inherited chromosomal imbalance":"ICB", "Digenic dominant":"DD"}
	url = "http://v1.marrvel.org/data/OMIM/"
	req = requests.get(url, params = {"geneSymbol": gene_id})
	data=req.json()
	mimNumber="NOOMIM"
	pheno=[]
	if data!=None:
		if "mimNumber" in data:
			mimNumber=str(data["mimNumber"])
		if "phenotypes" in data:
			for el in data["phenotypes"]:
				if el["phenotypeInheritance"] in inh:
					gg=inh[el["phenotypeInheritance"]]
				else:
					gg="ND"
				if "phenotypeMimNumber" in el and "phenotype" in el:
					pheno.append([str(el["phenotypeMimNumber"]), el["phenotype"], gg])
		if len(pheno)==0:
			pheno=["NOOMIM"]
	return mimNumber, pheno

def read_bd(bd, patern):
	"""read the animal ids and returns the entry of that id.
	if the id is not on the file, returns ND"""
	with open(bd) as f:
		line = next((l for l in f if patern in l), "ND")
	if len(line)==0:
		return "ND"
	else:
		return line


def make_it_shorter(l):
	gg=l.split(";")
	hh=len(gg)-1
	final=""
	while hh>0:
		if len(";".join(gg[:hh]))<=200:
			final=";".join(gg[:hh])+"..."
		hh-=1
	if final=="":
		return gg[0]+"..."
	else:
		return final

def get_data_animal(gene_id, bd):
	"""conects to the marrvel API and gets the animal models ID corresponding
	to the human gene ID. If it exists, returns the hyperlink ready for the
	output file, otherwise returns ND"""
	entrez=get_entrezid(gene_id)
	if entrez!="":
		bdd="http://api.marrvel.org/data/diopt/ortholog/gene/entrezId/"+str(entrez)
		req = requests.get(bdd)
		data=req.json()
		l=["flyBaseId","mgiId","wormBaseId", "zfinId", "rgdId"]
		final={"flyBaseId":"ND","mgiId":"ND","wormBaseId":"ND", "zfinId":"ND", "rgdId":"ND"}
		if len(data)>0:
			for el in data:
				if "gene2" in el and el['bestScore']==True:
					gg=el["gene2"]
					for database in l:
						if database in gg and final[database]=="ND":
							final[database]=str(gg[database])
			print("final_passo1", final)
			for key,value in final.items():
				if value=="ND":
					final[key]=["ND"]
				else:
					line=read_bd(bd[key][0],value)
					if line!="ND":
						rr=[]
						ll=line.split("\t")
						for el in ll[1].split(";"):
							if len(el)<=200:
								rr.append('=HYPERLINK("'+bd[key][1]+ll[0]+'","'+el+'")')
							else:
								tt=make_it_shorter(el)
								rr.append('=HYPERLINK("'+bd[key][1]+ll[0]+'","'+tt+'")')
						final[key]=rr
					else:
						final[key]=['=HYPERLINK("'+bd[key][1]+value+'","No disease association")']	
		else:
			final={"flyBaseId":["ND"],"mgiId":["ND"],"wormBaseId":["ND"], "zfinId":["ND"], "rgdId":["ND"]}
	else:
		final={"flyBaseId":["ND"],"mgiId":["ND"],"wormBaseId":["ND"], "zfinId":["ND"], "rgdId":["ND"]}
	return final


def get_entrezid(syn):
	for gene_id in syn:
		ii=gene_id.split("/")[0]
		bd="http://api.marrvel.org/data/gene/taxonId/9606/symbol/"+ii.strip()	
		print(bd)
		req= requests.get(bd)
		data=req.json()
		if 'entrezId' in data:
			return data['entrezId']
	else:
		print(gene_id)
		return ""

#final={"flyBaseId":["ND"],"mgiId":["ND"],"wormBaseId":["ND"], "zfinId":["ND"], "rgdId":["ND"]}
def get_all_animal_models(gene_id, params, syns):
	"""runs the get_data_animal function for all the animal models used.
	returns a list with all the results."""
	syns.append(gene_id)
	final=get_data_animal(syns, params["animals"])
	if final["wormBaseId"]==["ND"]:
		final["wormBaseId"]=['=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'+organism%3Aelegans&sort=score","'+gene_id+' C. elegans Uniprot entry")']
	else:
		final["wormBaseId"].append('=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'+organism%3Aelegans&sort=score","'+gene_id+' C. elegans Uniprot entry")')
	if final["flyBaseId"]==["ND"]:
		final["flyBaseId"]=['=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'+organism%3A%22Drosophila+melanogaster+%28Fruit+fly%29+%5B7227%5D%22&sort=score","'+gene_id+' Drosophila Uniprot entry")']
	else:
		final["flyBaseId"].append('=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'+organism%3Aelegans&sort=score","'+gene_id+' Drosophila Uniprot entry")')
	if final["mgiId"]==["ND"]:
		final["mgiId"]=['=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'&fil=organism%3A%22Mus+musculus+%28Mouse%29+%5B10090%5D%22&sort=score","'+gene_id+' Mouse Uniprot entry")']
	else:
		final["mgiId"]==["ND"].append('=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'&fil=organism%3A%22Mus+musculus+%28Mouse%29+%5B10090%5D%22&sort=score","'+gene_id+' Mouse Uniprot entry")')
	if final["rgdId"]==["ND"]:
		final["rgdId"]=['=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'&fil=organism%3A%22Rattus+norvegicus+%28Rat%29+%5B10116%5D%22&sort=score","'+gene_id+' Rat Uniprot entry")']
	else:
		final["rgdId"].append('=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'&fil=organism%3A%22Rattus+norvegicus+%28Rat%29+%5B10116%5D%22&sort=score","'+gene_id+' Rat Uniprot entry")')
	if final["zfinId"]==["ND"]:
		final["zfinId"]=['=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'&fil=organism%3A%22Danio+rerio+%28Zebrafish%29+%28Brachydanio+rerio%29+%5B7955%5D%22&sort=score","'+gene_id+' Zebrafish Uniprot entry")']
	else:
		final["zfinId"].append('=HYPERLINK("https://www.uniprot.org/uniprot/?query='+gene_id+'&fil=organism%3A%22Danio+rerio+%28Zebrafish%29+%28Brachydanio+rerio%29+%5B7955%5D%22&sort=score","'+gene_id+' Zebrafish Uniprot entry")')
	return [final["wormBaseId"], final["flyBaseId"], final["mgiId"], final["rgdId"],final["zfinId"]] 

	
def read_gtex_loops(gene_id, bd, t, chrr):
	"""reads the gtex/loops file and returns a list with all the entries
	of that respective gene. for gtex the gene must be in ensemblID and for
	theloops should be in gene name."""
	lines=["ND"]
	with open(bd) as f:
		if t=="gtex":
	  		lines = [l.split("\t")[-1].strip() for l in f if gene_id in l]
		else:
		 	lines = [l.split("\t")[-1].strip() for l in f if ("\t"+gene_id+" " in l and "chr"+chrr+"\t" in l)]
	return lines

def make_biblio_link_DD(ttype,band):
	return ('=HYPERLINK("https://pubmed.ncbi.nlm.nih.gov/?term=('+ttype+')+AND+('+band+')&sort=&filter=hum_ani.humans","'+band+" "+ttype+'")')

def make_biblio_link(geneid, li):
	"""makes the pubmed biblio link to be put on the output file"""
	gg="("+geneid+")"
	for el in li:
		if el!="nan":
			gg+="+OR+("+el+")"
	return ('=HYPERLINK("https://pubmed.ncbi.nlm.nih.gov/?term='+gg+'&sort=&filter=hum_ani.humans","'+geneid+'")')

