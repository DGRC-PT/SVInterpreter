#!/usr/bin/python3
encoding="UTF-8"
from collections import OrderedDict
from decimal import Decimal
import cProfile, pstats, io

def select_all_snps(infile, chrr, start, end, ont_file):
	"""Le o ficheiro da lista dos snps, e a regiao de interesse
	e retorna um dicionario com respectivo numero
	de SNPs como key, e como value uma lista de listas, em que
	cada lista corresponde a um trait e um p-value.
	Este metodo deve ser utilizado quando ja temos as regioes 
	intergenes definidas
	dic[#snps]=[traitA [p-valueA], traitB [p-valueB]]"""
	f=open(infile)
	dic={}
	for i in f:
		line=i.split("\t")
		if line[0]=="chr"+chrr or line[0]==chrr:
			if int(start)<=int(line[1]) and int(end)>=int(line[2]):
				if float(line[7].strip())<=0.00000005:
					if change_trait_format(line[3]) not in dic:
						dic[change_trait_format(line[3])]=[1,[float(line[7].strip())], 0, []]#primeiro os mais significativos, segudo os menos
					else:
						dic[change_trait_format(line[3])][0]+=1
						dic[change_trait_format(line[3])][1].append(float(line[7].strip()))#confirmar que isto funciona, sendo que o p-value esta em notacao cientifica
				elif float(line[7].strip())<0.00001:
					if change_trait_format(line[3]) not in dic:
						dic[change_trait_format(line[3])]=[0, [], 1,[float(line[7].strip())]]
					else:
						dic[change_trait_format(line[3])][2]+=1
						dic[change_trait_format(line[3])][3].append(float(line[7].strip()))#confirmar que isto funciona, sendo que o p-value esta em notacao cientifica					
	f.close()
	newdic={}
	for key,value in dic.items():
		if len(value[1])!=0:
			aa=formatnote(str(min(value[1])))
		if len(value[3])!=0:
			bb=formatnote(str(min(value[3])))
		if value[0]!=0 and str(value[0])+"b" not in newdic:
			newdic[str(value[0])+"b"]=[[key],["["+aa+"]"]]#confirmar que isto funciona, sendo que o p-value esta em notacao cientifica
		elif value[0]!=0 and str(value[0])+"b" in newdic:
			newdic[str(value[0])+"b"][0].append(key)
			newdic[str(value[0])+"b"][1].append("["+aa+"]")
		if value[2]!=0 and str(value[2]) not in newdic:
			newdic[str(value[2])]=[[key],["["+bb+"]"]]#confirmar que isto funciona, sendo que o p-value esta em notacao cientifica
		elif value[2]!=0 and str(value[2]) in newdic:
			newdic[str(value[2])][0].append(key)
			newdic[str(value[2])][1].append("["+bb+"]")
	ffdic=check_corres(ont_file, newdic)
	return ffdic


def formatnote(vv):
	if "E" not in vv:
		gg='%.2E' % Decimal(vv)
	else:
		gg=vv
	if "-0" in gg:
		aa=gg.replace("-0","-")
	else:
		aa=gg
	if ".00" in gg:
		bb=aa.replace(".00","")
	else:
		bb=aa
	return bb


def morestuf(l):
	n=[]
	for el in l:
		for r in el.split(","):
			n.append(r)
	return filter(None, n)

def get_correspondance(infile):
	f=open(infile)
	hpo={}#["",doid,efo,mesh]
	doid={}
	efo={}
	for i in f:
		aa=i.strip()
		line=aa.split("\t")
		if line[0]!="":
			for ele in line[0].split(","):
				hpo[ele.strip()]=morestuf(line[1:])
		elif line[1]!="":
			for ele in line[1].split(","):
				doid[ele.strip()]=morestuf(line[2:])
		elif line[2]!="" and line[3]!="":
			for ele in line[2].split(","):
				efo[ele.strip()]=line[3].split(",")
	f.close()
	return hpo, doid, efo

def desenvolve(val):
	n=[]
	for el in val:
		n.append(el.split("_")[0])
	return n

def check_corres(ont_corr, newdic):
	hpo, doid, efo=get_correspondance(ont_corr)
	fdic={}
	for key,value in newdic.items():
		fdic[key]=value
		l=[]
		for el in value[0]:
			if el.split("_")[0] in hpo:
				hh=desenvolve(value[0])
				for t in hpo[el.split("_")[0]]:
					if t in hh:
						aa=hh.index(t)
						l.append(aa)
			if el.split("_")[0] in doid:
				hh=desenvolve(value[0])
				for t in doid[el.split("_")[0]]:
					if t in hh:
						aa=hh.index(t)
						l.append(aa)
			if el.split("_")[0] in efo:
				hh=desenvolve(value[0])
				for t in efo[el.split("_")[0]]:
					if t in hh:
						aa=hh.index(t)
						l.append(aa)
		if len(l)>0:
			l.sort(reverse=True)
			for xx in l:
				del fdic[key][0][xx]
				del fdic[key][1][xx]
	ffdic=format_the_dict(fdic)
	return ffdic			

def format_the_dict(fdic):
	ffdic={}
	for key,value in fdic.items():
		aa=0
		ss=""
		while aa<len(value[0]):
			if ss=="":
				ss=value[0][aa]+" "+value[1][aa]
			else:
				ss+="; "+value[0][aa]+" "+value[1][aa]
			aa+=1
		ffdic[key]=ss
	return ffdic


def change_trait_format(trait):#only for the snps table
	"""Transforma o formato que vem por defeito na base de dados
	para a designacao dos traits, para a nossa forma:
	Transforma:
	DOID_3312_bipolar_disorder
	em
	DOID:3312_Bipolar_disorder
	Entra string sai string """
	dd=trait.split("_",2)
	iid=dd[0]+":"+dd[1]
	name=dd[2].capitalize()
	gg=name.replace("_"," ")
	return iid+"_"+gg


def get_cord(coo):
	co=coo.strip(")")
	c=co.split(":")
	coord=c[1].split("-")
	return [int(coord[0]), int(coord[1])]

def check_tad_disrupt(line, start, end):
	coord =get_cord(line[2])
	snps=0
	if (coord[0]<int(start) and coord[1]>int(start) and "(" in line[-1]) or (coord[0]<int(end) and coord[1]>int(end) and "(" in line[-1]):
		for el in line[-1].split(";"):
			if el!="\n":
				c=get_cord(el.split("(")[1])
				if c[0]>int(start) and c[0]<int(end):
					snps+=1
	elif int(line[6].strip())==0:
		snps="1"
	else:
		snps=line[6].strip()
	return str(snps)

def check_extra_tad(infile, ont_file, ensemblid, tadcord, is_tadstart):
	f=open(infile)
	f.readline()
	f.readline()
	dic={}
	for i in f:
		line=i.split("\t")
		if ensemblid.split(".")[0]==line[1].split(".")[0]:
			snps=0
			if line[6].strip()!="0":
				for el in line[-1].split(";"):
					if el!="\n":
						c=get_cord(el.split("(")[1])
						if is_tadstart==True and c[0]<int(tadcord):
							snps+=1
						if is_tadstart==False and c[0]> int(tadcord):
							snps+=1
			if float(line[5])<=0.00000005 and snps!=0:
				gg=str(snps)+"b"
			elif float(line[5])<0.00001 and snps!=0:
				gg=str(snps)
			elif snps==0:
				gg=""
			if gg!="":
				if gg.strip() not in dic:
					dic[gg.strip()]=[[line[3]+"_"+line[4].capitalize()],["["+line[5]+"]"]]
				else:
					dic[gg.strip()][0].append(line[3]+"_"+line[4].capitalize())
					dic[gg.strip()][1].append("["+line[5]+"]")
	f.close()
	ffdic=check_corres(ont_file, dic)
	nsn=[]
	if len(ffdic)!=0:
		nsn=list(ffdic.keys())
	return ffdic, nsn			
						
	

def read_gene_traits(infile, ensemblid,ont_file, start_tad,end_tad):
	"""Le o ficheiro que relaciona os genes com os snps, e retorna retorna 
	um dicionario com respectivo numero
	de SNPs como key (relativos ao gene, de acordo com o ensembl id,
	e como value uma lista de listas, em que
	cada lista corresponde a um trait e um p-value.
	dic[#snps]=[[traitA, p-valueA], [traitB, p-valueB]]"""
	f=open(infile)
	f.readline()
	f.readline()
	dic={}
	gg=""
	for i in f:
		line=i.split("\t")
		if ensemblid.split(".")[0]==line[1].split(".")[0]:
			snps=check_tad_disrupt(line, start_tad, end_tad)
			if float(line[5])<=0.00000005 and snps!="0":
				gg=snps+"b"
			elif float(line[5])<0.00001 and snps!="0":
				gg=snps
			elif snps=="0":
				gg=""
			if gg!="":
				if gg.strip() not in dic:
					dic[gg.strip()]=[[line[3]+"_"+line[4].capitalize()],["["+line[5]+"]"]]
				else:
					dic[gg.strip()][0].append(line[3]+"_"+line[4].capitalize())
					dic[gg.strip()][1].append("["+line[5]+"]")
	f.close()
	ffdic=check_corres(ont_file, dic)
	nsn=[]
	if len(ffdic)!=0:
		nsn=list(ffdic.keys())
	return ffdic, nsn

def establishcord(vv):
	dd=vv.split(":")
	gg=dd[1].split("-")
	if len(gg)>1:
		return dd[0],int(gg[0].replace(",","")),int(gg[1].replace(",",""))
	else:
		return dd[0],int(gg[0].replace(",","")),0

def get_index_for_bp(indexx, TAD):
	aa=[]
	bb=[]
	for el in indexx:
		if el.startswith("Break"):
			aa.append(dd)
			dd=indexx.index(el)
		else:
			dd=indexx.index(el)
	return aa

def read_regions(tad, dic):
	newb=OrderedDict()
	is_first=True
	final_tad=True
	for key, value in dic.items():
		chrr, start,end=establishcord(value[0])
		if start<int(tad[0]):
			newb[key+",start"]=start
			newb[key+",end"]=end
		elif start>int(tad[0]) and is_first==True:
			is_first=False
			newb["tad,start"]=int(tad[0])
			newb[key+",start"]=start
			if key.startswith("Brea")==False:
				newb[key+",end"]=end
		elif end>int(tad[1]):
			final_tad=False
			newb[key+",start"]=start
			if key.startswith("Brea")==False:
				newb[key+",end"]=end
		else:
			newb[key+",start"]=start
			if key.startswith("Brea")==False:
				newb[key+",end"]=end
		is_first=False
	if final_tad==True:
		newb["tad,end"]=int(tad[1])
	return newb

def make_regions(newb, dic, chrr):
	aa=0
	bb=OrderedDict()
	kefs=list(newb.keys())
	while aa<len(kefs):
		name=kefs[aa].split(",")
		if name[0] in dic and name[1]=="start" and name[0].startswith("Brea")==False:
			bb[name[0]]=dic[name[0]]
		else:
			if name[0].startswith("Brea"):
				 bb[name[0]]=dic[name[0]]
			if aa+1<len(kefs):
				s=newb[kefs[aa]]
				e=newb[kefs[aa+1]]
				if e>s:
					bb["chr"+chrr+":"+str(s)+"-"+str(e)]=["chr"+chrr, s, e]
		aa+=1
	return bb





	
	
