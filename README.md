# HardWiredGenome
Gene - Gene interaction network across the human Genome


This is a Wiring circuit for the genome constructed from existing literature and databases forming a Wiring Diagram across the human genome. 




### Databases
We use data from the following databases. The format of each file and method to link across each is described below.


1. STRING - Is a protein - protein interaction network.
2. HURI - HURI is a Protein - Protein interaction database 
3. HI-UNION - all HuRI interactions combined with CCSB screening efforts
4. LIT-BM - Records of High quality PPI from literature
5. Human TF - A list of all Human Transcription Factors from literature

#### Database file formats

```
data
└───STRING
│   │   protein.aliases.v12.0.txt
│   │   protein.links.v12.0.txt 
│   │   protein.info.v12.0.txt 
└───HURI
│   │   HuRI.tsv
│   │   HI-union.tsv
│   │   Lit-BM.tsv   
└───HTF
│   │   DatabaseExtract_v_1.01.csv
```

<b>STRING</b>
- STRING/protein.aliases.v12.0.txt - 3889207 lines + header - *string_protein_id* *alias* *source* 
- STRING/protein.links.v12.0.txt - 13715404 lines + header - *Protein1* *Protein2* *combined_score* - Contains the interactions between proteins 
- STRING/protein.info.v12.0.txt -  19699 lines + header - *string_protein_id*,*preferred_name* ... - Contains a list of all the proteins within the string database

<b>HuRI</b>
- HURI/HuRI.tsv - 52548 lines - *Protein1* *Protein2* - Has two proteins in each line denoting an interaction 
- HURI/HI-union.tsv - 64006 lines - *Protein1* *Protein2* -  Has two proteins in each line denoting an interaction
- HURI/Lit-BM.tsv - 13441 lines - *Protein1* *Protein2* - Has two proteins in each line denoting an interaction 


<b>Human TF</b>
2766 lines
- HTF/DatabaseExtract_v_1.01.csv  - 2765 lines + header - Ensembl ID, Is TF?
