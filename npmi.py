from __future__ import division
from Bio import Entrez
from Bio import Medline
import sys
import numpy as np
import time
import datetime

##########################################################################################
################################ Usage: python npmi.py ###################################
########################### developed by Gourab Das ######################################
###### For more inforrmation please visit https://github.com/dasgoba/npmistroke/ #########
############## Any query related to script please write to gourabdas0727@gmail.com #######
##########################################################################################
##########################################################################################

############# nPMI calculation script ###########

Entrez.email = "Email"	      ###### Please provide valid Email
Entrez.tool = "Script Name"   ###### Please provide script name



def allCount():
 strall="all[sb] AND hasabstract[text]"
 handle = Entrez.esearch(db="pubmed", term=strall) #, mindate=1811, maxdate=2010)
 record = Entrez.read(handle, validate=False)
 count_all = int(record["Count"]) ##### N= Total records in pubmed
 return(count_all)



def dfPubAbs(strp, count_all):
 stgr = strp+ "[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]"
 handle = Entrez.esearch(db="pubmed", term=stgr) #, mindate=1811, maxdate=2010)
 record = Entrez.read(handle, validate=False)
 count1 = int(record["Count"])
 p_gene = count1/count_all
 return(count1, p_gene)



def dfMutPubAbs(strp,count_all):

######### List of queries: Selected query should be uncommented before running the script ########

 ########### stgr = "(stroke[TIAB] OR Cerebrovascular[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]"    ########### Query Parent
 ########### stgr = '("Intracerebral hemorrhage"[TIAB] OR "Hemorrhagic Stroke"[TIAB] OR "Subarchanoid hemorrhage"[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child
 ########### stgr = '("Ischemic Stroke"[TIAB] OR (TOAST[TIAB] AND (Classification[TIAB] OR Subtypes[TIAB]))) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child
 ########### stgr = '"Large artery atherosclerosis"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Small vessel disease"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Cardioembolic disease"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Other determined etiology"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Undetermined etiology"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Intracerebral hemorrhage"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child
 ########### stgr = '"Subarachnoid Hemorrhage"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child

 strc = strp+"[TIAB] AND "+stgr
 handle = Entrez.esearch(db="pubmed", term=strc) #, mindate=1811, maxdate=2010)
 record = Entrez.read(handle, validate=False)
 count2 = int(record["Count"])
 pmu_gene = count2/count_all
 return(count2, pmu_gene)



def dfAnoPubAbs(count_all):

######### List of queries: Selected query should be uncommented before running the script ########

 ########### stgr = "(stroke[TIAB] OR Cerebrovascular[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]"     ########### Query Parent
 ########### stgr = '("Intracerebral hemorrhage"[TIAB] OR "Hemorrhagic Stroke"[TIAB] OR "Subarchanoid hemorrhage"[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child
 ########### stgr = '("Ischemic Stroke"[TIAB] OR (TOAST[TIAB] AND (Classification[TIAB] OR Subtypes[TIAB]))) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'    ########### Query Child
 ########### stgr = '"Large artery atherosclerosis"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Small vessel disease"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Cardioembolic disease"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Other determined etiology"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Undetermined etiology"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Grand Child
 ########### stgr = '"Intracerebral hemorrhage"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child
 ########### stgr = '"Subarachnoid Hemorrhage"[TIAB] AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]'   ########### Query Child

 handle = Entrez.esearch(db="pubmed", term=stgr) #, mindate=1811, maxdate=2010)
 record = Entrez.read(handle, validate=False)
 count3 = int(record["Count"])
 ano_gene = count3/count_all
 return(count3, ano_gene)



def nPMI(p1,p2,p12):
 pmi = np.ma.log(p12/(p1*p2)) ### Natural Logarithm and masked i.e. omits zero entries
 npmi = pmi/(-np.ma.log(p12)) #### Normalization done Bouma paper, 2007
 return(npmi)


#################### Calling Functions ################

now = datetime.datetime.now()
print ("Start date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))


count_all = allCount()                  #### count_all
y, ano_gene = dfAnoPubAbs(count_all)	#### p(y)

###### Please provide input and output files with paths
with open("Input File", 'r') as f, open("Output File", 'a') as fout:  ###### Kindly provide the inputs

    
 for strp in f:
  x, p_gene = dfPubAbs(strp,count_all)	#### p(x)
  xy, pmu_gene = dfMutPubAbs(strp,count_all)  #### p(x),p(y)
  #print p_gene,pmu_gene,ano_gene
  if(p_gene != 0):
   npmi = nPMI(p_gene,ano_gene,pmu_gene)
####   fout.write(strp.rstrip()+"	"+str(npmi))
   fout.write(strp.rstrip()+"	"+str(count_all)+"	"+str(x)+"	"+str(y)+"	"+str(xy)+"	"+str(npmi))
   fout.write("\n")
  else:
#### fout.write(str(strp.rstrip()+"	"+"NA"))
   fout.write(strp.rstrip()+"	"+str(count_all)+"	"+str(x)+"	"+str(y)+"	"+str(xy)+"	"+str(npmi))

   fout.write("\n")

  time.sleep(1)

now = datetime.datetime.now()
print ("End date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))
