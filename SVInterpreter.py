#!/usr/bin/python3
import cgi
import os
import cgitb
from time import time, localtime, strftime
import datetime
import calendar

cgitb.enable()

clock=strftime("%a, %b %d %Y %H:%M:%S", localtime())

def index():
	""" Show the default page
	"""
	print ('content-type: text/html')
	#print ('</html>')

index()


def showForm1():
	"""Show a form
	"""
	print ("""
<html>
    <head>
        <link rel="stylesheet" href="/test.css">
        <script src="/JS_scriptV2.js"></script>
        <link href="/nouislider.min.css" rel="stylesheet">
    </head>
    <div class="container">
        <form id="regForm" method=POST action="svinterpreter_aux/results_fantomV2.py">
            <script src="/nouislider.min.js"></script>
            <h1>Structural Variant Interpreter - SVInterpreter</h1>
            <a4>This tool was developed to support prediction of the phenotypic outcome of chromosomal or genomic structural variants (unbalanced and balanced translocations, inversion, insertion, deletions or duplications).</a4>
            <p></p>
            <a3> Please fill the following form with all the information about the structural variant to be analysed and respective phenotypic characteristics (optional). A table with relevant information for the evaluation of the structural variant will be retrived. </a3>
            <div class="row">
                <h4><center>Reference Human Genome (version)</h4>
                <div class="input-group">
				    <select name="version">
					    <option value="B" class="default" >Select Genome Version</option>
						<option value="hg19" ><b>Hg19</b></option>
						<option value="hg38" ><b>Hg38</b></option>
					</select>
                </div>
            </div>
            <div class="row">
			    <div class="input-group">
                    <h4><center>Cell line Hi-C data to use as reference</h4>
                    <a3><center>This data will be used to define the Topological Associated domains (TADs) boundaries and chromatin loops.</center></a3>
                    <a3><center>All data was retrived from <a href="http://3dgenome.fsm.northwestern.edu/publications.html">YUE Lab website.</a></center></a3>
                    <p></p>
                    <select name="tad">
                        <option value="B" class="default" >Select Cell-line</option>
                        <option value="consensus">Consensus TADs (Lifei 2019)</option>
                        <option value="IMR90">IMR90 (Rao 2014)</option>
                        <option value="LCL">LCL (Rao 2014)</option>
                        <option value="hesc">hESC (Dixon 2015)</option>
                        <option value="a549">A549 (Encode 2016)</option>
                        <option value="aorta">Aorta (Leung 2015)</option>
                        <option value="cortex">Cortex (Schmitt 2016)</option>
                        <option value="blader">Bladder (Schmitt 2016))</option>
                        <option value="lung">Lung (Schmitt 2016)</option>
                        <option value="huvec">HUVEC (Rao 2014)</option>
                        <option value="k562">K562 (Rao 2014)</option>
                    </select>
                </div>
            </div>
            <div class="row"> 
                <p><h4><center> Phenotypic description using HPO (optional)</h4></p>
                <div class="input-group">
                    <a3><center>The terms are separated by commas.</center></a3>
                    <p></p>
                    <input type="text" placeholder="HP:0000202, HP:0000157, HP:0006483, HP:0001640, HP:0001961,..." name="hpo_des"/>
                </div>
            </div>
            <div class="row">
                <h4><center> Highlighted Inheritance (optional)</h4></center>
                <a3><center>All phenotypes are analyzed and presented, but only the ones with the user-selected inheritance are highlighted on the output.</center></a3>
                <p></p>
                <div class="input-group">
				    <select name="inh">
				        <option value="B" class="default" >Select Inheritance</option>
					    <option value="AD">Autossomal Dominant (AD)</option>
						<option value="AR">Autossomal Recessive (AR)</option>
						<option value="PD">Pseudoautosomal Dominant (PD)</option>
						<option value="PR">Pseudoautosomal Recessive (PR)</option>
						<option value="DD">Digenic Dominant (DD)</option>
						<option value="DR">Digenic Recessive(DR)</option>
						<option value="IC">Isolated Cases (IC)</option>
						<option value="ICB">Inherited chromosomal imbalance (ICB)</option>
						<option value="Mu">Multifactorial(Mu)</option>
						<option value="SMo">Somatic mosaicism (SMo)</option>
						<option value="SMu">Somatic mutation (SMu)</option>
						<option value="XL">X-linked (XL)</option>
						<option value="XLD">X-linked Dominant (XLD)</option>
						<option value="XLR">X-linked Recessive (XLR)</option>
						<option value="YL">Y-linked (YL)</option>
					</select>
                </div>
            </div>
            <div class="row">
                <div class="input-group">
				    <h4><center> Type of structural variant</h4>
                    <select id="seltest", name="tt", onchange="yesnoCheck(tt);">
                        <option value="B" class="default" >Select Structural variant</option>
                        <option value="Balanced_translocation">Balanced Translocation</option>
                        <option value="Unbalanced_translocation">Unbalanced Translocation</option>
                        <option value="Inversion">Inversion</option>
                        <option value="Deletion">Deletion</option>
                        <option value="Duplication">Duplication</option>
                        <option value="Insertion">Insertion</option>
                        <option value="Spec_Rg">Query Genomic region</option>
                    </select>
                </div>
            </div>
            <div class="row">
                <div class="ifYes", id="ifYes">
                    <a3>This tools accepts coordinates or intervals, with or without commas.</a3>
                    <p></p>
                    <a4>Chromosome A</a4>
                    <input id="chrA" type="text" class="bal" placeholder="6" name="chrA"/>
                    <p></p>
                    <a4>Chromosome B</a4>
                    <input id="chrB" type="text" class="bal" placeholder="7" name="chrB"/>
                    <p></p>
                    <a4>Breakpoint A</a4>
                    <input id="brA" type="text" class="bal" placeholder="168,529,498" name="brA"/>
                    <p></p>
                    <a4>Breakpoint B</a4>
                    <input id="brB" type="text" class="bal" placeholder="168,529,498" name="brB"/>
                </div>
            </div>
            <div class="row">
                <div class="ifYes", id="ifno">
                    <a3>This tools accepts coordinates or intervals, with or without commas.</a3>
                    <p></p>
                    <a4>Chromosome </a4>
                    <input id="chrA" type="text" class="bal" placeholder="6" name="chrA"/>
                    <p></p>
                    <a4>Region (start-end)</a4>
                    <input id="brA" type="text" class="bal" placeholder="168,529,498-168,540,877" name="brA"/>
                </div>
            </div>
            <div class="row">
                <div class="ifYes", id="ifdd">
                    <a3>This tools accepts coordinates or intervals, with or without commas.</a3>
                    <p></p>
                    <a4>Chromosome </a4>
                    <input id="chrA" type="text" class="bal" placeholder="6" name="chrA"/>
                    <p></p>
                    <a4>Region (start-end)</a4>
                    <input id="brA" type="text" class="bal" placeholder="168,529,498-168,540,877" name="brA"/>
                    <p></p>
                    <input class="checkbox_input" type="checkbox" name="coupon_question" value="1" id="isov">
                    <label for="isov">Overlap with public Databases?</label>
                </div>
            </div>
            <div class="row">
                <div class="ifYes", id="typeofanaly">
				    <h4><center> Analysis based on</h4></center>
				    <select name="tttype" onchange="shohide(this)">
				        <option value="B" class="default" >Select strategy </option>
					    <option value="ggg">TADs</option>
						<option value="spec">A specific region</option>
					</select>
                </div>
            </div>
            <div class="row">
                <div class="ifYes" id="formTADs">
					<h4><center> TADs to analyse</h4></center>
					<a3><center>Use the sliders to define the TADs to analyse. By default, only the breakpoint TAD (brTAD) is analysed.</center></a3>
				    <p></p>
                    <script src="nouislider.min.js"></script>
                    <div id="slider">
                    </div>
                    <input type="hidden" name="slider_control" id="slider_control" value="" />
                    <p>&nbsp;</p>
                    <script>
                        var slider = document.getElementById('slider');
                        noUiSlider.create(slider, {
                            start: [-0.3, 0.3],
                            connect: true,
                            range: {
                                'min': -5,
                                'max': 5
                            },
                            pips: {
                                mode: 'positions',
                                stepped: true,
                                density: 10,
                                values: [0,10,20,30,40,50,60,70,80,90,100],
                                format: {
                                    to: function(value) {
                                        if (value<0){
                                            return ["TAD-1","TAD-2","TAD-3", "TAD-4", "TAD-5"][(value*-1)-1];
                                        }else{
                                            return ["brTAD", "TAD+1", "TAD+2", "TAD+3", "TAD+4", "TAD+5"][value];
                                        }
                                    },
                                    from: Number
                                }
                            }
                        });
                        var slider_control = document.getElementById('slider_control');
                        slider.noUiSlider.on('update', function( values) {
                            slider_control.value = values.join(' - ');
                        });
                    </script>
                </div>
            </div>
            <div class="row">
                <div class="ifYes" id="dvPassport">
                    <h4>Insert the region to analyse</h4>
                    <a3>The Region must follow the chr:start-end format.</a3>
                    <a3>Make sure to insert a region that includes the variant breakpoint.</a3>
                    <p></p>
                    <input type="text" placeholder="chr5:23345543-23357777" name="specreg"/>
                </div>
            </div>
            <div class="row">
                <div class="ifYes" id="dvPassport2">
                    <h4>Insert the regions to analyse</h4>
                    <a3>The Regions must follow the chr:start-end format.</a3>
                    <a3>Insert the regions corresponding to each chromossome of the rearrangement.</a3>
					<a3>Make sure to insert regions that include the variant breakpoints.</a3>
                    <p></p>
					<a4>Chromosome A region </a4>
                    <input type="text" placeholder="chr5:23345543-23357777" name="specreg1"/>
					<a4>Chromosome B region </a4>
					<input type="text" placeholder="chr7:44345543-55777766" name="specreg2"/>
                </div>
            </div>
            <div class="row">
                <div  style="display: none;" id="ifoverlap">
                    <h4><center>Databases to used on the overlap search</center></h4>
                    <input type="checkbox" name="dats[]" value="DGV" id="dgv" checked />
                    <label style="word-wrap:break-word" for="dgv">DGV database</label>
                    <p><input type="checkbox" name="dats[]" value="1000Genomes" id="1000genomes" checked/>
                    <label style="word-wrap:break-word" for="1000genomes">1000 Genomes database</label>
                    <p><input type="checkbox" name="dats[]" value="ClinGen" id="clingen" checked/>
                    <label style="word-wrap:break-word" for="clingen">ClinGen/ClinVar databases</label>
                    <p><input type="checkbox" name="dats[]" value="deldupsindrome" id="deldupsindrome" checked/>
                    <label style="word-wrap:break-word" for="deldupsindrome">Deletion/Duplication syndromes</label>
                    <p><input type="checkbox" name="dats[]" value="CoeCoop" id="coecoop" checked/>
                    <label style="word-wrap:break-word" for="coecoop">DDelay studies (Cooper et al. 2011; Coe et al. 2012)</label>
                    <p><input type="checkbox" name="dats[]" value="collins" id="collins" checked/>
                    <label style="word-wrap:break-word" for="collins">Collins et al. 2017</label>
                    <p><input type="checkbox" name="dats[]" value="chaisson" id="chaisson" checked/>
                    <label style="word-wrap:break-word" for="chaisson">Chaisson et al. 2019</label>
                    <p><input type="checkbox" name="dats[]" value="gnomad" id="gnomad" checked/>
                    <label style="word-wrap:break-word" for="gnomad">Gnomad SV database</label>
                    <p>&nbsp;</p>
                    <h4><center>Type of overlap search:</center></h4>
                    <select id="ovl" name="ovl">
                        <option value="B" class="default" >Selection strategy</option>
                        <option value="mutual">Mutual Overlap</option>
                        <option value="full">Query comprised by the reference</option>
                    </select>
                </div>
            </div>
            <div class="row">
                <div  style="display: none;" id="ifmutual">
                    <h4> <center> Mutual Overlap </center></h4>
                    <a3>Choosing Mutual Overlap, this tool will apply a mutual overlap cut-off on the search:</a3>
                    <a3> Database hits are only retrieved if percentages of overlap of query vs database and database vs query are above the cut-off.</a3>
                    <a3><p><b>Overlap cutoff (1-100)%:</b></p></a3>
                    <input type="text" name="perc" value="70">
                </div>
            </div>
            <div class="row">
                <div  style="display: none;" id="iffull">
                    <h4> <center> Query comprised by the reference </center></h4>
                    <a3> Choosing query comprised by the reference, this tool will retrieve all database hits that covers 100% of the query, independently of the database entry size.<a3>   
                </div>
            </div>
            <div class="row">
                <div class="ifYes", id="ifins">
                    <a3>This tools accepts coordinates or intervals, with or without commas.</a3>
                    <p></p>
                    <a4>Recipient Chromosome </a4>
                    <input id="chrA" type="text" class="bal" placeholder="6" name="chrA"/>
                    <p></p>
                    <a4>Donor Chromosome </a4>
                    <input id="chrB" type="text" class="bal" placeholder="8" name="chrB"/>
                    <p></p>
                    <a4>Recipient Breakpoint</a4>
                    <input id="brA" type="text" class="bal" placeholder="116812107-116912603" name="brA"/>
                    <p></p>
                    <a4>Inserted region</a4>
                    <input id="brB" type="text" class="bal" placeholder="168,539,498" name="brB"/>
                </div>
            </div>
            <div class="row">
                <div class="ifYes", id="ifmeh">
                    <a3>This tools accepts coordinates or intervals, with or without commas.</a3>
                    <p></p>
                    <a4>Chromosome A </a4>
                    <input id="chrA" type="text" class="bal" placeholder="6" name="chrA"/>
                    <p></p>
                    <a4>Breakpoint A</a4>
                    <input id="brA" type="text" class="bal" placeholder="116812107-116912603" name="brA"/>
                    <p></p>
                    <input type="radio" name="vvv" value="del" id="del1"/>
                    <label for="del1">Deletion</label>
                    <input type="radio" name="vvv" value="dup" id="dup1"/>
                    <label for="dup1">Duplication</label>
                    <p></p>
                    <a4>Chromosome B </a4>
                    <input id="chrB" type="text" class="bal" placeholder="8" name="chrB"/>
                    <p></p>
                    <a4>Breakpoint B</a4>
                    <input id="brB" type="text" class="bal" placeholder="16812107-16912603" name="brB"/>
                    <p></p>
                    <input type="radio" name="ddd" value="del" id="del"/>
                    <label for="del">Deletion</label>
                    <input type="radio" name="ddd" value="dup" id="dup"/>
                    <label for="dup">Duplication</label>
                </div>
            </div>
            <p><input type="submit" value="Submit"></b></p></center>
        </form>
        <p><button onClick="window.location.href=window.location.href">Clean Form</button></p>
    </div>
<p></p>
<a3><center>If you using this tool please acknowledge it by citing <a href="https://link.springer.com/article/10.1007/s00439-020-02121-x">our reference publication</a></center>
<center><address>
Correspondance: <a href="mailto:doencasgenomicas@insa.min-saude.pt">Genomic Diseases Group</a>
</address>
<center><aaa><a href="http://www.insa.min-saude.pt/category/areas-de-atuacao/genetica-humana/">Department of Human Genetics</a></aaa></center>
<p>National Institute of Health Doutor Ricardo Jorge</p> </aaa></center>
<center><img src="https://cld.pt/dl/download/bf231ea4-336c-47c2-98a9-5129c3af3510/aaa.png" width=500 height=80 border=0 alt=""><br><br />
<center><p><rodape>This file was last modified 28/12/2020</p></font></a3>
</html>
""")

showForm1()
#<input type="reset">
#<p><button onClick="window.location.href=window.location.href">Clean Form</button></p>
