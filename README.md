# A text data driven method for stroke associated gene identification

Association between stroke related keywords and human gene symbols have been measured from the PubMed using normalized pointwise mutual information (nPMI). PubMed advanced queries have been prepared utilizing the stroke related keywords from domain knowldge and gene symbols. For each type and sub-type of the stroke and each gene symbol, PubMed title and abstracts have been searched with the query and no. of hits (document frequency DF) have been recorded. DFs have been converted to nPMI using the formula by Bouma et al. [Proceedings of GSCL 30 (2009): 31-40]. A python script has been written for nPMI calculation.

Prerequisites: 
        1. Biopython (>= 1.65)
        2. NumPy ()
        3. A text file of gene symbols. One symbol per line.
        4. A output text file
        5. PubMed advanced query
 
Compatibility:
        python 2 or 3

Usage: python npmi.py

Sample Inputs:

Filename: input.txt

IMPACT
MTHFR
APOE
TNF
ACE
SET
NOTCH3

Sample advanced query:

Query for stroke: "(stroke[TIAB] OR Cerebrovascular[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]"

Sample Output:

Gene_symbol | #N | #X | #Y | #XY | nPMI
--- | --- | --- | --- |--- |--- 
IMPACT | 20809291 | 73226 | 9623 | 419 | 0.233135901


	     N	       X	     Y	    XY	     nPMI
	      20809291	 73226	  9623  	419	   0.233135901
MTHFR	        20816630	 4297	    9626	  307	   0.453021431
APOE	        20802743	 6763	    9619	  267    0.394728553
TNF	          20816630	 31754	  9626	  261	   0.255329041
ACE	          20802743	 5515   	9619	  238	   0.399504136
SET	          20816630	 75862	  9626	  192	   0.145395645
NOTCH3	      20816630	 859	    9626	  191	   0.532774728

where,

N: Total no. of PubMed entires on that date
X: No. of hits related to query 'x'
Y: No. of hits related to query 'y'
XY: No. of hits related to query X combined with query Y
nPMI: 




