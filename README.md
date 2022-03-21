## A novel text data driven method for stroke associated gene identification

Association between stroke related keywords and human gene symbols have been measured from the PubMed using normalized pointwise mutual information (nPMI). PubMed advanced queries have been prepared utilizing the stroke related keywords from domain knowldge and gene symbols. For each type and sub-type of the stroke and each gene symbol, PubMed title and abstracts have been searched with the query and no. of hits (document frequency DF) have been recorded. DFs have been converted to nPMI using the formula by Bouma et al. [Proceedings of GSCL 30 (2009): 31-40]. A python script has been written for nPMI calculation.

<br/>

**Prerequisites:** <br/>
 
1. Biopython (>= 1.65) <br/>
2. NumPy (latest) <br/>
3. A text file of gene symbols. One symbol per line <br/>
4. A output text file <br/>
5. PubMed advanced query <br/>
 
<br/>

**Compatibility:**
        python 2 or 3

<br/>

**Usage:  python npmi.py**

<br/>

**Sample Inputs:**
Filename: input.txt
--- | 
IMPACT 
MTHFR 

<br/>

**Sample advanced query:**

1. Query for stroke: <br/>
"(stroke[TIAB] OR Cerebrovascular[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]" 
2. Query for Hemorrhagic Stroke:<br/>
'("Intracerebral hemorrhage"[TIAB] OR "Hemorrhagic Stroke"[TIAB] OR "Subarchanoid hemorrhage"[TIAB]) AND (gene[TIAB] OR genes[TIAB]) AND hasabstract[text]' 


<br/>

**How to run the script:**

1. Open npmi.py in suitable text editor (preferably Notepad++ or EditPlus). <br/> 
2. Enter email id and script name in the specified positions (Follow comments in the script). <br/>
3. Uncomment/ Write required PubMed advanced query in the specified positions (Follow comments in the script. See samples above or in the script). <br/>
4. Please enter input and output text file with absolute paths in the specified positions (Follow comments in the script). <br/>
5. Save the modified script. <br/>
6. Run by typingthe following command in linux terminal or windows dos promt: **python npmi.py** <br/>

<br/>


**Sample Output:**

Gene_symbol | #N | #X | #Y | #XY | nPMI
--- | --- | --- | --- |--- |--- 
IMPACT | 20809291 | 73226 | 9623 | 419 | 0.233135901
MTHFR |	20816630 | 4297 | 9626 | 307 | 0.453021431


where,

N: Total no. of PubMed entires on that date <br/>
X: No. of hits related to query 'x' <br/>
Y: No. of hits related to query 'y' <br/>
XY: No. of hits related to query X combined with query Y <br/>
nPMI: Normalized pointwise mutual information calculated using Bouma et al.

<br/>

**Please Cite:** <br/>
Gourab, Das, and Pradeep Gupta. "Potential Key Genes Associated with Stroke types and its subtypes: A Computational Approach." Neuroscience Informatics (2022): (Submitted).

<br/>

**Troubleshooting:**

This script calculates nPMI using NCBI The Entrez Programming Utilities (E-utilities). Hence, kindly prepare the input file with at max 10 gene symbols at a time for quick output and avoid PubMed database connetion error. If still encounter PubMed connetion error, please re-run the script. For more details please follow NCBI https://www.ncbi.nlm.nih.gov/books/NBK25497/. 

<br/>

**Contact:**

Dr. Gourab Das <br/>
Email: gourabdas0727@gmail.com

