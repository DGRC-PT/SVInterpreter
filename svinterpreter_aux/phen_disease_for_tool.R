#input is a case phenotype and a OMIM/ORPHA disease.
##################################################################################
#INPUT
rm(list=ls())
.libPaths( c( .libPaths(), "/home/dgrc/R/x86_64-pc-linux-gnu-library/3.6") )
args <- commandArgs(TRUE)
#indicate the path to the rda file
load(file="hpodata.rda")
#casephen<-strsplit(args[1], ",")
#casephen<-c("HP:0011344", "HP:0000750", "HP:0001344", "HP:0002136", "HP:0001290", "HP:0011968", "HP:0001250","HP:0012443", "HP:0000708","HP:0000729","HP:0040082","HP:0040196","HP:0000280","HP:0430028","HP:0005280","HP:0000463","HP:0000453","HP:0000154","HP:0010800","HP:0010804","HP:0000215","HP:0000687","HP:0008897","HP:0000486","HP:0030680","HP:0001629","HP:0001636","HP:0012020","HP:0000032","HP:0001741","HP:0000034","HP:0001537","HP:0002558","HP:0002943","HP:0001776")
#entrez geneID to be compared with the case phenotype
geneID<-args[length(args)]
#geneID<-"OMIM:617616"
#number of random observations for calculating the p-value. The more the better, but it will increase the time span
iterations<-100
#p-value cutoff
pvalue<-0.05
##################################################################################
casephen<-vector()
aa=1
while (aa<length(args)){
  casephen <-append(casephen,as.character(args[aa]))
  aa=aa+1
}
###################################################################
print(casephen)
print(geneID)
#Run HPOSim
library(HPO.db)
library(HPOSim)
.initialize()
combinemethod = "funSimMax"
method = "Resnik"
IC <- get("termIC", envir = HPOSimEnv)
#input the case phenotype

#calculate max score possible
Terms1<-casephen
Terms2<-casephen
Terms1 <- RemoveTermsWithoutIC(Terms1, "PA", IC)
Terms2 <- RemoveTermsWithoutIC(Terms2, "PA", IC)
Max_Ker <- getTermListSim(Terms1, Terms2, combinemethod, method, IC, verbose)

#make the similarity between the gene of interest and the case phenotype
Terms2<-as.character(unlist(disease[geneID]))
print(geneID)
print(Terms2)
if (length(Terms2)!=0){
  Terms1 <- RemoveTermsWithoutIC(Terms1, "PA", IC)
  Terms2 <- RemoveTermsWithoutIC(Terms2, "PA", IC)
  Case_Ker <- getTermListSim(Terms1, Terms2, combinemethod, method, IC, verbose)
  
  #to modulate the p-value, we have to overlap the case phenotype with all possible phenotypes available
  v<-vector()
  bb=1
  while (bb<=iterations){
    Terms1 <- sample(hpo_terms, length(casephen))
    Terms2<-as.character(unlist(disease[geneID]))
    Terms1 <- RemoveTermsWithoutIC(Terms1, "PA", IC)
    Terms2 <- RemoveTermsWithoutIC(Terms2, "PA", IC)
    Ker <- getTermListSim(Terms1, Terms2, combinemethod, method, IC, verbose)
    v<-append(v,Ker)
    bb=bb+1
  }
  
  qq<-length(v[v >= Case_Ker])
  if (qq==0){
    qq=qq+1
  }
  result<-(qq+1)/iterations
  print(paste("PhenSSc", as.character(signif(Case_Ker, digits=3))," (P=", as.character(signif(result, digits=3)), "; MaxSSc",as.character(signif(Max_Ker, digits=3)), ")"))
}else{
  print("Sc ND ; P = ND")
}

