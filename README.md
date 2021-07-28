#SVInterpreter

##This tool was developed to support prediction of the phenotypic outcome of chromosomal or genomic structural variants (unbalanced and balanced translocations, inversion, insertion, deletions or duplications). And is available online in: [Here](https://dgrctools.insa.min-saude.pt/)

###It is possible to run SVInterpreter locally using a preexistant Apache / cgi-bin configuration.

##Pre-Requisites
* Apache configuration
* Python collections
* Python openpyxl
* Python sys
* Python subprocess
* Python cgi
* Python time
* Python random
* Python pickle
* Python re
* Python pybiomart
* Python pandas
* Python cProfile
* Python pstats
* Python io
* Python requests
* Python urllib
* Python json
* Python decimal
* Python bisect
* R
* R httr package
* R hpo.bd and HPOSim packages. For more information see [Here](https://github.com/DGRC-PT/HPOSim_Helper_Scripts)

##Resources on Data directory
###This distribution of SVInterpreter comes with several data, including:
* Chromosome sizes for Hg19 and Hg38
* Cytoband for [Hg19](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1134961021_TAIMkBe3Bbnq2IN3vGYza4zjOVB0&db=hg19&c=chr1&g=cytoBand) and [Hg38](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1134961021_TAIMkBe3Bbnq2IN3vGYza4zjOVB0&db=hg38&c=chr1&g=cytoBand)
* Centromere for [Hg19](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1134961021_TAIMkBe3Bbnq2IN3vGYza4zjOVB0&db=hg19&g=gap) and [Hg38](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1134961021_TAIMkBe3Bbnq2IN3vGYza4zjOVB0&db=hg38&c=chr1&g=centromeres)
* [ACMG reportable genes](https://www.coriell.org/1/NIGMS/Collections/ACMG-59-Genes)
* GeneHancer cluster of interactions for [Hg19](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1134961949_gPA9W50onLggaQbKpU4x2AnnvVPP&db=hg19&g=geneHancer) and [Hg38](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1134961949_gPA9W50onLggaQbKpU4x2AnnvVPP&db=hg38&g=geneHancer)
* [TADs and Loops for Hg19 and Hg38](http://3dgenome.fsm.northwestern.edu/publications.html)
* [Haploinsuficiency](https://www.deciphergenomics.org/about/downloads/data) and [Triplosensitivity](https://dosage.clinicalgenome.org/help.shtml)
* [Infertility data](https://pubmed.ncbi.nlm.nih.gov/30865283/)
* [Observed vs expected Lof](https://gnomad.broadinstitute.org/)
* [PannelApp data](https://panelapp.genomicsengland.co.uk/)

###Since some resources used by SVInterpreter come for diferent databases that should be regularly updated, we show below how to compile the data needed by SVInterpreter.
* ###Flybase data
    1 - Go to [Flybase](http://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release)and download the "disease_model_annotation" file

    2 - Decompress the file

    3 - Use the script prepare_flybase.py to create the bd file
 `python3 prepare_flybase.py disease_model_anotation > flybase_formated_bd`

    4 - Place the output file inside the "data" directory

* ###MGI data

    1 - Go to [MGI](http://www.informatics.jax.org/downloads/reports/index.html#pheno) and download the "MGI_GenePheno.rpt" and "VOC_MammalianPhenotype.rpt" files

    2 - Decompress the files if needed

    3 - Use the script prepare_flybase.py to create the bd file
 `python prepare_mgi.py MGI_GenePheno.rpt VOC_MammalianPhenotype.rpt > mgi_formated_bd`

    4 - Place the output file inside the "data" directory

* ###RGD data

    1 - Go to [RGD](https://download.rgd.mcw.edu/data_release/) and download the "GENES_RAT" file

    2 - Decompress the files if needed

    3 - Use the script prepare_rat.py to create the bd file
 `python prepare_rat.py GENES_RAT > rat_formated_bd`

    4 - Place the output file inside the "data" directory

* ###Wormbase data

    1 - Go to [wormbase](https://wormbase.org/about/userguide#3--10) and download the "disease_association.WS276.daf.txt" file

    2 - Download the [DOID.obo](https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo) file

    3 - Use the script prepare_rat.py to create the bd file
 `python prepare_wormbase.py disease_association.WS276.daf.txt doid.obo > wormbase_formated_bd`

    4 - Place the output file inside the "data" directory


* ###ZFin data

    1 - Go to [Zfin](https://zfin.org/downloads) and download the "Gene to disease via Ontology" file

    2 - Use the script prepare_rat.py to create the bd file
 `python prepare_zfin.py gene2DiseaseViaOntology > zfin_formated_bd`

    4 - Place the output file inside the "data" directory

* ###ClinGen data
You can use the data from [Here](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=phenDis&hgta_track=clinGenComp&hgta_table=clinGenGeneDisease&hgta_doSchema=describe+table+schema) and create a tab-separated text file with the following columns: Gene symbol, Disease name, MONDO_ID, HGNC_ID, Classification, link to the clinGen database. Place the file in the  "data" directory

* ####DDG2P Data
Download the DDG2P data from [Here](https://www.deciphergenomics.org/about/downloads/data), decompress the file, and place it on the "data" directory.


###Instalation
Simply paste the SVInterpreter.py and svinterpreter_aux folder on your cgi-bin folder.


###Notes:
* It is possible that some paths defined on the scripts need to be edited to work in your specific workspace
