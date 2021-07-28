#!/usr/bin/env python3
import sys
import subprocess
import cgi
import cgitb;cgitb.enable()
import time
import random
import string
import comand_fantomV2 as comand_fantom
import subprocess
###HTML dealli

def make_table(outfile1, chrA, chrB, brA, brB, tt, der, link, name_3):
    aa=open("/var/www/html/dgrctools.insa.min-saude.pt/public_html/outputs/out_debug", "w")
    aa.write(";".join(['python3','report_table_fantom.py', outfile1, chrA, chrB, brA, brB, tt, der, link, name_3]))
    #aa.close()
    #p = subprocess.Popen(['python3','report_table_fantomV2.py', outfile1, chrA, chrB, brA, brB, tt, der, link, name_3], stdout=subprocess.PIPE, stderr=subprocess.STDOUT); print('   ')
    p = subprocess.Popen(['python3','report_table_fantomV2.py', outfile1, chrA, chrB, brA, brB, tt, der, link, name_3], stdout=aa, stderr=aa); print('   ')
    aa.close()

def gg(comand):
    print('content-type:text/html; charset=utf-8 \n\n')
    print('<html><head><title>HTML Meta Tag</title><meta http-equiv = "refresh" content = "1; url = '+comand+'"></head><body></body></html>')

#aa=open("/var/www/html/dgrctools.insa.min-saude.pt/public_html/outputs/out_debug", "w")

cgitb.enable()
form=cgi.FieldStorage()
l=''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
link="/var/www/html/dgrctools.insa.min-saude.pt/public_html/outputs/out_"+l+".html"
comand_fantom.index2(link)
params, outfile1, chrA, chrB, brA, brB, tt, der, name_3=comand_fantom.conr(form, link)
#aa.write(chrA)
if chrA!="aa":
	comand_fantom.remake(link)
	comand_fantom.save_obj(params, link.split("/")[-1])
	gg("/outputs/out_"+l+".html")
	#aa.close()
	make_table(outfile1, chrA, chrB, brA, brB, tt, der, link, name_3)
else:
	gg("/outputs/out_"+l+".html")

#aa.close()
