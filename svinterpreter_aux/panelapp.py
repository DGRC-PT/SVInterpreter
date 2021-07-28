#!/usr/bin/python3
encoding="UTF-8"
from sys import argv
import urllib.request
import urllib.error
import json
import string
from collections import OrderedDict

def get_genes_panel(panel):
	request =("https://panelapp.genomicsengland.co.uk/WebServices/get_panel/"+panel)
	try:
		response = urllib.request.urlopen(request)
		kittens = response.read()
		bb=json.loads(kittens.decode('utf-8'))
		aa=bb["result"]
		for ii in aa["Genes"]:
			for key,value in ii.items():
				if value==None:
					ii[key]="nd"
			print(panel.replace("%20"," ")+"\t"+ii["GeneSymbol"]+"\t"+ii["ModeOfInheritance"]+"\t"+ii["Penetrance"]+"\t"+ii["LevelOfConfidence"])
	except urllib.error.HTTPError as e:
		print (e)


def run_panels(infile):
	f=open(infile)
	for i in f:
		gg=i.strip()
		hh=gg.replace(" ","%20")
		get_genes_panel(hh)
	f.close()


#run_panels(argv[1])


def try_api(gene):
	path=OrderedDict({"HighEvidence":[],"ModerateEvidence":[], "LowEvidence":[]})
	request =("https://panelapp.genomicsengland.co.uk/WebServices/search_genes/"+gene+"/")
	try:
		response = urllib.request.urlopen(request)
		kittens = response.read()
		bb=json.loads(kittens.decode('utf-8'))
		for ek in bb["results"]:
			if 'LevelOfConfidence' in ek:
				if ek['LevelOfConfidence'] not in path:
					path[ek['LevelOfConfidence']]=[ek['SpecificDiseaseName']]
				else:
					path[ek['LevelOfConfidence']].append(ek['SpecificDiseaseName'])
		return path
	except urllib.error.HTTPError as e:
		return (e)



def testlistgenes(infile):
	f=open(infile)
	for i in f:
		try_api(i.strip())
	f.close()

#testlistgenes(argv[1])


def try_panel(panel):
        request =("https://panelapp.genomicsengland.co.uk/WebServices/get_panel/"+panel)
        try:
                response = urllib.request.urlopen(request)
                kittens = response.read()
                bb=json.loads(kittens.decode('utf-8'))
        except urllib.error.HTTPError as e:
                print (e)


def get_panel(panel):
	#request =("https://panelapp.genomicsengland.co.uk/WebServices/list_panels/?Name="+panel)
	request =("https://panelapp.genomicsengland.co.uk/WebServices/get_panel/"+panel)
	try:
		response = urllib.request.urlopen(request)
		kittens = response.read()
		bb=json.loads(kittens.decode('utf-8'))
		if "result" not in bb:
			return("None")
		else:
			ek=bb["result"]
			if ek["DiseaseGroup"]!="":
				return(ek["DiseaseGroup"])
			else:
				return("None")
	except urllib.error.HTTPError as e:
		return (e)

def testlistpanel(infile):
	f=open(infile)
	for i in f:
		aa=i.strip()
		get_panel(aa.replace(" ","%20"))
	f.close()


def manipulate_genes(gene):
	cats={"Actionable information":"A", "Cancer Programme":"B", "Cardiovascular disorders":"C", "Ciliopathies":"D", 
	"Dermatological disorders":"E", "Dysmorphic and congenital abnormality syndromes":"F", "Endocrine disorders":"G",
	"Gastroenterological disorders": "H", "Growth disorders":"I", "Haematological and immunological disorders":"J",
	"Haematological disorders":"K", "Hearing and ear disorders":"L", "Metabolic disorders":"M", "Neurology and neurodevelopmental disorders":"N",
	"Ophthalmological disorders":"O", "Renal and urinary tract disorders": "P", "Respiratory disorders":"Q", "Rheumatological disorders":"R", 
	"Skeletal disorders":"S", "Tumour syndromes":"T", "None":"U"}
	result=try_api(gene)
	if len(result)!=0:
		out={"***":[], "**":[], "*":[]}
		for key,value in result.items():
			for ele in value:
				panel=get_panel(ele.replace(" ","%20"))
				if panel in cats:
					leter=cats[panel]
				else:
					for el in list(string.ascii_uppercase):
						if el not in cats.values():
							cats[panel]=el
							break
					leter=cats[panel]
				if key=="HighEvidence" and leter not in out["***"]:
					out["***"].append(leter)
				elif key=="ModerateEvidence" and leter not in out["***"] and leter not in out["**"]:
					out["**"].append(leter)
				elif key=="LowEvidence" and leter not in out["***"] and leter not in out["**"] and leter not in out["*"]:
					out["*"].append(leter)
		l=""
		if len(out["***"])!=0:
			out["***"].sort()
			l+="***"+",".join(out["***"])
		if len(out["**"])!=0:
			out["**"].sort()
			if len(l)>0:
				l+=";"
			l+="**"+",".join(out["**"])
		if len(out["*"])!=0:
			out["*"].sort()
			if len(l)>0:
				l+=";"
			l+="*"+",".join(out["*"])
		return(l)
	else:
		return("")

def read_and_transform(merged):
	f=open(merged, encoding='utf-8')
	dic={}
	for i in f:
		line=i.split("\t")
		entry1=line[0]+" / "+line[1]+"_"+line[2]
		entry11=line[4]+":"+line[3]+", "+line[-1].strip()+" cases"
		entry2=line[1]+" / "+line[0]+"_"+line[2]
		entry22=line[4]+":"+line[3]+", "+line[-1].strip()+" cases"
		if line[0] not in dic:
			dic[line[0]]=[['=HYPERLINK("'+line[5]+'","'+entry1+'")', '=HYPERLINK("'+line[5]+'","'+entry11+'")']]
		else:
			dic[line[0]].append(['=HYPERLINK("'+line[5]+'","'+entry1+'")', '=HYPERLINK("'+line[5]+'","'+entry11+'")'])
		if line[1] not in dic:
			dic[line[1]]=[['=HYPERLINK("'+line[5]+'","'+entry2+'")', '=HYPERLINK("'+line[5]+'","'+entry22+'")']]
		else:
			dic[line[1]].append(['=HYPERLINK("'+line[5]+'","'+entry2+'")', '=HYPERLINK("'+line[5]+'","'+entry22+'")'])
	f.close()
	return dic

def tiny_transform(v):
	f=open(v, encoding='utf-8')
	dic={}
	for i in f:
		line=i.split("\t")
		dic[line[0]]=[['=HYPERLINK("'+line[-1].strip()+'","'+line[1]+'")', '=HYPERLINK("'+line[-1].strip()+'","'+line[2]+'")']]
	f.close()
	return dic
