from sys import argv

from collections import OrderedDict


rr=OrderedDict([('chr1:72300527-72346279', [['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv24502&ref=GRCh37/hg19","esv24502")', '45687', 'arr[GRCh38]1p31.1(72300599-72346286)x3', 'gain+loss', '28/40 70.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/19812545","Conrad_et_al_2009")', '', '99.843', '99.985'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv33765&ref=GRCh37/hg19","esv33765")', '45147', 'arr[GRCh38]1p31.1(72300811-72345958)x3', 'gain+loss', '29/51 56.863%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/17666407","de_Smith_et_al_2007")', '', '98.678', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv3584698&ref=GRCh37/hg19","esv3584698")', '42212', 'arr[GRCh38]1p31.1(72303253-72345465)x3', 'gain', '10/34 29.412%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/24956385","Mokhtar_et_al_2014")', '', '92.263', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv1012282&ref=GRCh37/hg19","nsv1012282")', '59480', 'arr[GRCh38]1p31.1(72302803-72362283)x3', 'gain', '1/29084 0.003%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/25217958","Coe_et_al_2014")', '', '95.025', '73.093'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv10250&ref=GRCh37/hg19","nsv10250")', '46560', 'arr[GRCh38]1p31.1(72300436-72346996)x3', 'gain+loss', '24/31 77.419%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/18304495","Perry_et_al_2008")', '', '100.0', '98.265'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv2422081&ref=GRCh37/hg19","esv2422081")', '45063', 'arr[GRCh38]1p31.1(72300402-72345465)x1', 'deletion', '1027/1184 86.74%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/20811451","Altshuler_et_al_2010")', '', '98.221', '99.723'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv24502&ref=GRCh37/hg19","esv24502")', '45687', 'arr[GRCh38]1p31.1(72300599-72346286)x3', 'gain+loss', '28/40 70.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/19812545","Conrad_et_al_2009")', '', '99.843', '99.985'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv2480676&ref=GRCh37/hg19","esv2480676")', '46606', 'seq[GRCh38]1p31.1(72300197-72346803)x1', 'deletion', '1/1 100.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/19546169","McKernan_et_al_2009")', '', '100.0', '98.168'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv2579314&ref=GRCh37/hg19","esv2579314")', '46317', 'seq[GRCh38]1p31.1(72300059-72346376)x1', 'loss', '1/1 100.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/19546169","McKernan_et_al_2009")', '', '100.0', '98.78'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv2749808&ref=GRCh37/hg19","esv2749808")', '45553', 'seq[GRCh38]1p31.1(72300599-72346152)x1', 'deletion', '93/96 96.875%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/23290073","Wong_et_al_2012b")', '', '99.565', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv33765&ref=GRCh37/hg19","esv33765")', '45147', 'arr[GRCh38]1p31.1(72300811-72345958)x3', 'gain+loss', '29/51 56.863%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/17666407","de_Smith_et_al_2007")', '', '98.678', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv3563904&ref=GRCh37/hg19","esv3563904")', '52999', 'seq[GRCh38]1p31.1(72297318-72350317)x1', 'deletion', '1/767 0.13%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/23714750","Boomsma_et_al_2014")', '', '100.0', '86.326'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv4289&ref=GRCh37/hg19","esv4289")', '45759', 'seq[GRCh38]1p31.1(72300579-72346338)x1', 'loss', '1/1 100.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/18987735","Wang_et_al_2008")', '', '99.886', '99.871'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv9190&ref=GRCh37/hg19","esv9190")', '46030', 'seq[GRCh38]1p31.1(72300395-72346425)x1', 'loss', '0/1 0.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/19470904","Ahn_et_al_2009")', '', '100.0', '99.396'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=esv995847&ref=GRCh37/hg19","esv995847")', '46380', 'seq[GRCh38]1p31.1(72300094-72346474)x1', 'deletion', '1/3 33.333%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/20482838","Pang_et_al_2010")', '', '100.0', '98.646'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv10250&ref=GRCh37/hg19","nsv10250")', '46560', 'arr[GRCh38]1p31.1(72300436-72346996)x3', 'gain+loss', '24/31 77.419%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/18304495","Perry_et_al_2008")', '', '100.0', '98.265'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv437860&ref=GRCh37/hg19","nsv437860")', '40799', 'arr[GRCh38]1p31.1(72302069-72342868)x1', 'loss', '46/269 17.1%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/16468122","McCarroll_et_al_2006")', '', '89.174', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv442908&ref=GRCh37/hg19","nsv442908")', '42232', 'arr[GRCh38]1p31.1(72303233-72345465)x1', 'loss', '238/270 88.148%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/18776908","McCarroll_et_al_2008")', '', '92.306', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv498672&ref=GRCh37/hg19","nsv498672")', '45517', 'seq[GRCh38]1p31.1(72300641-72346158)x1', 'loss', '1/9 11.111%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/21111241","Kidd_et_al_2010b")', '', '99.486', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv508282&ref=GRCh37/hg19","nsv508282")', '56506', 'seq[GRCh38]1p31.1(72297409-72353915)x1', 'deletion', '1/4 25.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/20534489","Teague_et_al_2010")', '', '100.0', '80.968'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv511699&ref=GRCh37/hg19","nsv511699")', '45564', 'seq[GRCh38]1p31.1(72300609-72346173)x1', 'loss', '1/1 100.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/21212237","Arlt_et_al_2011")', '', '99.589', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv513990&ref=GRCh37/hg19","nsv513990")', '44264', 'arr[GRCh38]1p31.1(72300945-72345209)x1', 'loss', '416/2366 17.582%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/21397061","Campbell_et_al_2011")', '', '96.748', '100.0'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv546532&ref=GRCh37/hg19","nsv546532")', '50821', 'arr[GRCh38]1p31.1(72292938-72343759)x1', 'loss', '1/17421 0.006%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/21841781","Cooper_et_al_2011")', '', '94.492', '85.067'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv546549&ref=GRCh37/hg19","nsv546549")', '32740', 'arr[GRCh38]1p31.1(72300402-72333142)x1', 'loss', '5/17421 0.029%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/21841781","Cooper_et_al_2011")', '', '71.287', '99.618'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv819206&ref=GRCh37/hg19","nsv819206")', '45800', 'seq[GRCh38]1p31.1(72300532-72346332)x1', 'loss', '1/2 50.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/19587683","Kim_et_al_2009")', '', '99.989', '99.884'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv821382&ref=GRCh37/hg19","nsv821382")', '45874', 'seq[GRCh38]1p31.1(72300544-72346418)x1', 'deletion', '1/1 100.0%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/20802225","Ju_et_al_2010")', '', '99.963', '99.697'], ['chr1:72300527-72346279', '45752', '=HYPERLINK("http://dgv.tcag.ca/dgv/app/variant?id=nsv947431&ref=GRCh37/hg19","nsv947431")', '56152', 'seq[GRCh38]1p31.1(72287838-72343990)x1', 'deletion', '8/97 8.247%', '=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/23825009","Sudmant_et_al_2013")', '', '94.997', '77.402']])]) 

def av(v1,v2):
	aa=0.0
	aa=(float(v1)+float(v2))/2
	return aa
	
def perform_best_hit( rr):
	mean_overlap=OrderedDict()
	freq=OrderedDict()
	sample=OrderedDict()
	best_hit=OrderedDict()
	link=OrderedDict()
	l2={}
	to_pass=[]
	for key, value2 in rr.items():
		if len(value2)==0:
			best_hit[key]=""
			to_pass=""
		elif len(value2)==1:
			cc=value2[0][6].split(" ")
			fre=cc[1]
			nn=value2[0][2].split('"')[3]
			gg=value2[0][2].replace('"'+nn+'"', '"'+nn+' ('+fre+')"')
			best_hit[key]=gg
			to_pass.append(nn)#name
			to_pass.append(value2[0][2])#link
			to_pass.append(nn+" (freq. "+fre+"; %overlap "+value2[0][-1]+")")

		else:
			for value in value2:
				nn=value[2].split('"')[3]
				mean_overlap[nn]=av(value[-1], value[-2])
				cc=value[6].split(" ")
				fre=cc[1].replace("%","")
				freq[nn]=float(fre)
				sample[nn]=int(cc[0].split("/")[1])
				link[nn]=value[2].replace('"'+nn+'"', '"'+nn+' ('+fre+'%)"')
				l2[nn]=[nn, value[2], nn+" (freq. "+fre+"; %overlap "+value[-1]+")"]
			name, v=search_best_hit(mean_overlap, freq, sample, link)
			best_hit[key]=v
			to_pass=l2[name]
			mean_overlap=OrderedDict()
			freq=OrderedDict()
			sample=OrderedDict()
	return best_hit, to_pass
	
def write(best_hit):
	for key, value in best_hit.items():
		print(key+"\t"+value)

def get_val_from_val(dic, top):
	r=top*-1
	dd=[]
	indxs=set([])
	for key, value in dic.items():
		dd.append(value)
	gg=(sorted(dd)[r:])
	for el in gg:
		jj=gg.index(el)
		indxs.add(jj)
	return gg, indxs

def search_best_hit(mean_overlap, freq, sample,link):
	top=2
	while top<=len(mean_overlap):
		mo, moidx=get_val_from_val(mean_overlap, top)
		fo, foidx=get_val_from_val(freq, top)
		so, soidx=get_val_from_val(sample, top)
		finals=list(moidx&foidx&soidx)
		if len(finals)==1:
			name=list(mean_overlap.keys())[finals[0]]
			return(name, link[name])
			#return(name+" ("+str(freq[name])+"%)")
		if len(finals)>1 or top==len(mean_overlap):
			name=list(mean_overlap.keys())[finals[-1]]
			return(name, link[name])
			#return(name+" ("+str(freq[name])+"%)")
		if len(finals)==0:
			top+=1

#write(perform_best_hit(rr))
	
	
