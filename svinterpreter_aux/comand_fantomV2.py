#!/usr/bin/env python3
coding="UTF-8"

#report table subscript
##############################################################
###Dependencies:
import sys
import subprocess
from openpyxl import Workbook
from openpyxl.styles import Font, Fill, PatternFill
from openpyxl.styles.borders import Border, Side
import cgi
import cgitb;cgitb.enable()
from collections import OrderedDict
import time
import random
#import report_table_fantom as rt
import parameters_file
import cross_deletions_with_dgv_results_Vfor_new
###HTML dealling
import pickle
import re


def take_space_coma(val):## vai para o command.py
	aa=val.strip(" ")
	if ":" in aa:
		ss=aa.split(":")[1]
	else:
		ss=aa
	bb=ss.replace(",","")
	return bb

def take_space_coma_and_sep(val):## vai para o command.py
	aa=val.strip(" ")
	bb=aa.replace(",","")
	return bb.split("-")

def make_tads(v):
	"""takes the nouislider values in str (minval-maxval) and
	transforms it in ints in a list: [min_val, max_val]. If
	only one tad is selected, the list is retrived with only one element"""
	gg=v.rsplit("-",1)
	if round(float(gg[0]))!=round(float(gg[1])):
		return [round(float(gg[0])), round(float(gg[1]))]
	else:
		return [round(float(gg[0]))]

def make_title_tads(gg):
	"""gets the output of make_tads and transforms it into a text
	representation of the same. ex: input:[-2,4]; output:TAD-2 to TAD4.
	this information is used on the output page."""
	hh=[]
	for el in gg:
		if el==0:
			hh.append("brTAD")
		elif el>0:
			hh.append("TAD+"+str(el))
		else:
			hh.append("TAD"+str(el))
	return " to ".join(hh)

def make_reg(gg):
	dd=gg.split(":")
	cords=dd[1].split("-")
	return dd[0].replace("chr", ""), cords[0], cords[1]

def write_error(te, link):
	out=open(link, "a")
	out.write('<h3><center><p style="color:red"><b>Input error!!</h3></center></p></b>')
	out.write ('<a3><p><center><b>'+te+'</b></a3></p></center>')
	out.write('<input type="button" value="Click here to rectify your input and try again" onclick="history.back()">')
	out.write('</div>')
	out.write('''
<center>If you using this tool please acknowledge it by citing <a href="https://link.springer.com/article/10.1007/s00439-020-02121-x">our reference publication</a></center>
<center><address>
Correspondance: <a href="mailto:doencasgenomicas@insa.min-saude.pt">Genomic Diseases Group</a></center>
</address>
<center><aaa><a href="http://www.insa.min-saude.pt/category/areas-de-atuacao/genetica-humana/">Department of Human Genetics</a></aaa></center>
<center><aaa><p>National Institute of Health Doutor Ricardo Jorge</p> </aaa></center>
<center><img src="https://cld.pt/dl/download/bf231ea4-336c-47c2-98a9-5129c3af3510/aaa.png" width=500 height=80 border=0 alt=""><br><br /></center>
<center><p><rodape>This file was last modified 28/12/2020</p></font></center></html>''')
	out.close()


def validate_region(reg, ch, sz, tt):
	aa=True
	if re.search('[a-zA-Z]', reg)!=None:
		aa=False
	else:
		r=take_space_coma_and_sep(reg)
		if (tt=="Deletion" or tt=="Duplication" or tt=="Inversion" or tt=="Spec_Rg") and len(r)==1:
			aa=False
		if len(r)>1:
			if int(r[1])<int(r[0]):
				aa=False
			elif int(r[1])>sz[ch]:
				aa=False
		else:
			if int(r[0])>sz[ch]:
				aa=False
	return aa

def verify_input(form, link):
	final=True
	chrs=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19", "20", "21","22","X","Y"]
	if form["version"].value=="hg19":
		sz={"1":249250621, "2":243199373, "3":198022430, "4":191154276, "5":180915260, "6":171115067, "7":159138663, "8":146364022, "9":141213431, "X":155270560, "Y":59373566, "10":135534747, "11":135006516, "12":133851895, "13":115169878, "14":107349540, "15":102531392, "16":90354753, "17":81195210, "18":78077248, "19":59128983, "20":63025520, "21":48129895, "22":51304566}
	if form["version"].value=="hg38":
		sz={"1":248956422, "2":242193529, "3":198295559, "4":190214555, "5":181538259, "6":170805979, "7":159345973, "8":145138636, "9":138394717, "X":156040895, "Y":57227415, "10":133797422, "11":135086622, "12":133275309, "13":114364328, "14":107043718, "15":101991189, "16":90338345, "17":83257441, "18":80373285, "19":58617616, "20":64444167, "21":46709983, "22":50818468}
	if form["version"].value=="B":
		write_error("Please specify a genome version and try again.", link)
		final=False
	elif form["tad"].value=="B":
		write_error("Please specify a cell line to use as reference and try again.", link)
		final=False
	elif form["tt"].value=="B":
		write_error("Please specify a structual variant to analyse and try again.", link)
		final=False
	elif form["tttype"].value=="B":
		write_error("Please specify if the analysis should be based on TADs or a specific region, and try again.", link)
		final=False
	else:
		tt=form["tt"].value
		if  tt=="Deletion" or tt=="Duplication" or tt=="Inversion" or tt=="Spec_Rg":
			if form["chrA"].value not in chrs:
				write_error("Please correctly specify the chromosome (only the ID of the chromosome without any chr prefix, special characters or spaces) and try again.", link)
				final=False
			elif validate_region(form["brA"].value, form["chrA"].value, sz, tt)==False:
				write_error("Please specify a valid region to be analysed, and try again. Make sure that you did not used invalid characters, that your start coordinate < end coordinate, and that the region exists on the specified chromosome and genome version.", link)
				final=False
			elif "coupon_question" in form and form["coupon_question"].value=="1":
				if "dats[]" not in form:
					write_error("Please specify which databases to use on the overlap search, and try again.", link)
					final=False
				elif form["ovl"].value=="B":
					write_error("Please specify which databases search strategy should be used, and try again.", link)
					final=False
				elif form["ovl"].value=="mutual" and (int(form["perc"].value)<0 or int(form["perc"].value)>100):
					write_error("Please specify a valid overlap percentage cutoff (1-100) and try again.", link)
					final=False
				elif form["tttype"].value=="spec":
					chh, start,end=make_reg(form["specreg"].value)
					if chh!=form["chrA"].value :
						write_error("Please make sure that you specified specific regions to analyse that overlap the structural variant, and try again.", link)
						final=False
					if chh not in chrs or int(end)<int(start) or int(end)>sz[chh]:
						write_error("Please specify a valid specific region to be used as reference to the analysis, and try again.", link)
						final=False
		else:
			if form["chrA"].value not in chrs or form["chrB"].value not in chrs:
				write_error("Please correctly specify the chromosomes (only the ID of the chromosome without any chr prefix, special characters or spaces) and try again.", link)
				final=False
			elif validate_region(form["brA"].value, form["chrA"].value, sz, tt)==False or validate_region(form["brB"].value, form["chrB"].value, sz,tt)==False:
				write_error("Please specify a valid breakpoints, and try again. Make sure that you did not used invalid characters, that your start coordinate < end coordinate (in case of a region), and that the breakpoint exists on the specified chromosome and genome version.", link)
				final=False
			elif form["tttype"].value=="spec":
				if "specreg1" not in form or "specreg2" not in form:
					write_error("Please specify a valid specific regions to be used as reference to the analysis, and try again.", link)
					final=False
				else:
					chh, start,end=make_reg(form["specreg1"].value)
					chh2, start2,end2=make_reg(form["specreg2"].value)
					if chh!=form["chrA"].value or chh2!=form["chrB"].value:
						write_error("Please make sure that you specified specific regions to analyse that overlap the structural variant, and try again.", link)
						final=False
					elif chh not in chrs or int(end)<int(start) or int(end)>sz[chh] or chh2 not in chrs or int(end2)<int(start2) or int(end2)>sz[chh2]:
						write_error("Please specify a valid specific regions to be used as reference to the analysis, and try again.", link)
						final=False
	return final

def conr(form, link):
	result_ver=verify_input(form, link)
	if result_ver==False:
		return "aa", "aa", "aa", "aa", "aa", "aa", "aa", "aa", "aa"
	else:
		out=open(link, "a")
		alts=parameters_file.alts#carrega o dicionario com as alteracoes
		version=form["version"].value#c
		if version=="hg19":
			params=parameters_file.hg19
			al=parameters_file.hg19_lines
		else:
			params=parameters_file.hg38
			al=parameters_file.hg38_lines
		tdds=[]
		ta=[]
		params["hits"]=[]
		params["hpo_des"]=""
		if "hpo_des" in form:
			if len(form["hpo_des"].value)>1:
				gg=form["hpo_des"].value
				uu=""
				for el in gg.split(","):
					uu+=' "'+el.strip()+'"'
				params["hpo_des"]=uu
		if form["tt"].value=="Spec_Rg":
			thebps=take_space_coma_and_sep(form["brA"].value)
			params["tads"]=[[form["chrA"].value, thebps[0],thebps[1]]]
			title=''
			params["ttype"]="spec"
			t2='<a3><p><center><b>Analysis of: </b> Genomic Region </a3></p></center>'
		else:
			if form["tttype"].value=="ggg":
				params["tads"]=make_tads(form["slider_control"].value)#gets the tads for the params dic
				title='<a3><p><b><center>TADs to analyse: </b>'+ make_title_tads(params["tads"])+'</p></center></a3>'
			else:
				if "specreg" in form:
					chh, start,end=make_reg(form["specreg"].value)
					if form["tt"].value!="Inversion":
						params["tads"]=[[chh, start,end]]
					elif form["tt"].value=="Inversion":
						thebps=take_space_coma_and_sep(form["brA"].value)
						params["tads"]=[[chh, str(min(int(thebps[0]),int(start))), str(max(int(thebps[0]),int(start)))],[chh, str(min(int(thebps[1]),int(end))), str(max(int(thebps[1]),int(end)))]]
					title='<p><a3><center><b>Region to analyse: </b>'+ chh +':'+str(start)+'-'+str(end)+'</p></center></a3>'
				elif "specreg1" in form and "specreg2" in form:
					chh, start,end=make_reg(form["specreg1"].value)
					chh2, start2,end2=make_reg(form["specreg2"].value)
					params["tads"]=[[chh, start,end], [chh2, start2,end2]]	
					title='<a3><p><b><center>Regions to analyse: </b>'+ chh +':'+str(start)+'-'+str(end)+';'+chh2 +':'+str(start2)+'-'+str(end2)+'</p></center></a3>'
			t1=form["tt"].value
			t2='<a3><p><center><b>Type of alteration: </b>'+ t1.replace("_", " ")+'</a3></p></center>'
			params["ttype"]=form["tttype"].value
		tadd=form["tad"].value#linha
		params["tad"]=al[str(tadd)]
		tt=form["tt"].value
		params["tt"]=tt
		if "inh" in form:
			params["input_inh"]=form["inh"].value
		if "ddd" in form:
			if form["ddd"].value =="A":
				der=form["chrA"].value
			else:
				der=form["chrB"].value
		else:
			der=""
		out.write('<br /><a3><font size="5em;"><center><b> Input parameters: </center></b></font></a3>')
		out.write ('<a3><p><center><b>Genome version: </b>'+ version+'</a3></p></center>')
		out.write ('<a3><p><center><b>Reference: </b>'+ params["tad"][0]+'</a3></p></center>')
		if "hpo_des" not in form:
			hpo="None"
		else:
			hpo=form["hpo_des"].value
		out.write ('<a3><p><center><b>HPO Description: </b>'+hpo+'</a3></p></center>')
		if params["input_inh"]!="B":
			out.write ('<a3><p><center><b>Selected Inheritance: </b>'+params["input_inh"]+'</a3></p></center>')
		else:
			out.write ('<a3><p><center><b>Selected Inheritance: </b> None </a3></p></center>')
		out.write(title)
		out.write (t2)
		if  tt=="Deletion" or tt=="Duplication" or tt=="Inversion" or tt=="Spec_Rg":
			chh=form["chrA"].value 
			name11=time.strftime(alts[tt]+chh+"_"+params["version"]+"_"+tadd+"_%d-%m-%Y_report_table.xlsx")
			outfile1=("/var/www/html/dgrctools.insa.min-saude.pt/public_html/outputs/"+name11)
			out.write("<a3><p><center><b>Chromosome: </b>"+chh+"</a3></p></center>")
			out.write("<a3><p><center><b>Region: </b>"+form["brA"].value+"</a3></p></center>")
			if tt=="Inversion":
				aa=chh
				thebps=take_space_coma_and_sep(form["brA"].value)
				#name1, name2= rt.main(params, outfile1, chh, aa, thebps[0],thebps[1], tt, der)##########################
				out.close()
				return params, outfile1, chh, aa, thebps[0],thebps[1], tt, der, ""
			#elif tt=="Spec_Rg":
			#	thebps=take_space_coma_and_sep(form["brA"].value)
				#out.write("checgou aui")
				#out.write(tt)
			#	out.close()
			#	return params, outfile1, chh, "del", thebps[0],thebps[1], tt, der, ""
			else:
				aa="del"
				if "coupon_question" in form and form["coupon_question"].value=="1":
					dats=[]
					fd={"DGV":"DGV", "1000Genomes":"1000Genomes", "ClinGen":"ClinGen", "deldupsindrome":"Deletion and Duplication syndromes", "CoeCoop":"DDelay studies", "collins":"Collins et al. 2017", "chaisson":"Chaisson et al. 2019", "gnomad":"GnomadSV"}
					for i in form["dats[]"]:
						dats.append(fd[i.value])
					if form["ovl"].value=="mutual":
						ovl="Mutual overlap"
					else:
						ovl="Query comprised reference"
					out.write ('<a3><p><center><b>Databases overlap: </b> '+",".join(dats)+'</a3></p></center>')
					out.write ('<a3><p><center><b>Overlap Strategy: </b> '+ovl+'</a3></p></center>')
					if ovl=="Mutual overlap":
						perc=form["perc"].value
						out.write ('<a3><p><center><b>Overlap cutoff: </b> '+perc+'%</a3></p></center>')
					name_3, hits=cross_deletions_with_dgv_results_Vfor_new.exect(form)###############################
					out.write("Chegou aqui!")
					thebps=take_space_coma_and_sep(form["brA"].value)
					params["hits"]=hits###############################
					#name1, name2= rt.main(params, outfile1, chh, aa, thebps[0],thebps[1], tt, der)
					out.close()
					return params, outfile1, chh, aa, thebps[0],thebps[1], tt, der, name_3
					#print("<aaa><center><a href=/outputs/"+name_3+">Download CNV table!</a></center></aaa>")####################################
				else:
					out.write ('<a3><p><center><b>Databases overlap: </b> None </a3></p></center>')
					out.write ('<a3><p><center><b>Overlap Strategy: </b> None </a3></p></center>')
					thebps=take_space_coma_and_sep(form["brA"].value)
					#name1, name2= rt.main(params, outfile1, chh, aa, thebps[0],thebps[1], tt, der)############################
					out.close()
					return params, outfile1, chh, aa, thebps[0],thebps[1], tt, der, ""
		else:
			chrA=form["chrA"].value
			chrB=form["chrB"].value
			name11= time.strftime(alts[tt]+chrA+"_"+chrB+"_"+form["version"].value+"_"+tadd+"_%d-%m-%Y_report_table.xlsx")
			params["html_A"]=time.strftime(alts[tt]+chrA+"_%d-%m-%Y_.html")
			params["html_B"]=time.strftime(alts[tt]+chrB+"_%d-%m-%Y_.html")
			outfile1=("/var/www/html/dgrctools.insa.min-saude.pt/public_html/outputs/"+name11)
			out.write("<a3><p><center><b>Chromosome A: </b>"+chrA+"</a3></p></center>")
			out.write("<a3><p><center><b>Chromosome B: </b>"+chrB+"</a3></p></center>")
			out.write("<a3><p><center><b>Breakpoint A: </b>"+take_space_coma(form["brA"].value)+"</a3></p></center>")
			out.write("<a3><p><center><b>Breakpoint B: </b>"+take_space_coma(form["brB"].value)+"</a3></p></center>")#######################
			#name1, name2= rt.main(params, outfile1, chrA, chrB, take_space_coma(form["brA"].value), take_space_coma(form["brB"].value), tt, der)
			out.close()
			return params, outfile1, chrA, chrB, take_space_coma(form["brA"].value), take_space_coma(form["brB"].value), tt, der, ""
			#print ('<a3><font size="5em;"><center><p><b>Output:</center></b></p></a3></font>')
			#print("<a3><p><center><b>Rearrangement A: </b>"+name1+"</p></center></a3>")
			#print("<a3><p><center><b>Rearrangement B: </b>"+name2+"</p></center></a3>")
			#print('<center><p><a href=/'+params["html_A"]+name1+' Online vizualization</center></a></p>')
			#print('<center><p><a href=/'+params["html_B"]+name2+' Online vizualization</center></a></p>')
		#print('<center><p><a href=/outputs/'+name11+'>Download report table!</center></a></p>')

def save_obj(obj, g):
	name=g[:-5]
	with open("/var/www/html/dgrctools.insa.min-saude.pt/public_html/outputs/"+name + '.pkl', 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def remake(link):
	out=open(link,"a")
	out.write('<h3><center><p style="color:red"><b>SVInterpreter is creating your table. The running time depends on the complexity of the table and the size of the region. Please be patient....</h3></center></p></b>')
	out.write('</div>')
	out.write('''
<center>If you using this tool please acknowledge it by citing <a href="https://link.springer.com/article/10.1007/s00439-020-02121-x">our reference publication</a></center>
<center><address>
Correspondance: <a href="mailto:doencasgenomicas@insa.min-saude.pt">Genomic Diseases Group</a></center>
</address>
<center><aaa><a href="http://www.insa.min-saude.pt/category/areas-de-atuacao/genetica-humana/">Department of Human Genetics</a></aaa></center>
<center><aaa><p>National Institute of Health Doutor Ricardo Jorge</p> </aaa></center>
<center><img src="https://cld.pt/dl/download/bf231ea4-336c-47c2-98a9-5129c3af3510/aaa.png" width=500 height=80 border=0 alt=""><br><br /></center>
<center><p><rodape>This file was last modified 28/12/2020</p></font></center></html>''')
	out.close()

def index2(link):
	out=open(link,"w")
	out.write ('''
<html>
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<link href="https://fonts.googleapis.com/css?family=Arial" rel="stylesheet">
<head>
<meta charset="UTF-8">
''')
	out.write('<meta http-equiv="refresh" content="10">')
	out.write('''
<title>SVInterpreter</title>
<link rel="stylesheet" href="/outputs/test.css">
<script src="/outputs/JS_script.js"></script>
<link href="/outputs/nouislider.min.css" rel="stylesheet">
</head>
<body>
<div class="container">
<br><br/>
<center><h1>SVInterpreter - Search Results</h1></center>
<hgh><center><p><a3>The retrieved data from each breakpoint is compiled in a table that can be vizualized on the browser or downloaded in xlsx.</a3></p></center></hgh>
<center><hgh><button onclick="goBack()">New search</button></center></hgh>
<script>
function goBack() {
  window.location.replace("/cgi-bin/SVInterpreterV2.py");
}
</script>
</div>
<div class="container">
''')
	out.close()

