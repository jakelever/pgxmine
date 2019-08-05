# Supplementary Materials

This file contains details of specific methods for the PGxMine project.

## Drug lists

The list of drugs to extract and map is created by the [createDrugList.py](https://github.com/jakelever/pgxmine/blob/master/createDrugList.py) script. This uses [MeSH](https://www.nlm.nih.gov/databases/download/mesh.html), [DrugBank](https://www.drugbank.ca/releases/latest) and [PharmGKB](https://www.pharmgkb.org/downloads) drug lists to create a mapping file from MeSH ID to PharmGKB ID. It also provides a normalized name using the DrugBank name, and flags whether the drug is a cancer treatment (if 'Antineoplastic' appears in a DrugBank category). 

It removes drugs in categories: 'Elements','Adenine Nucleotides','Sweetening Agents','Salt Solutions','Supplements','Solvents','Electrolyte Solutions','Food Additives','Food','Lactates','Diluents','Gases','Mineral Supplements','Basic Lotions and Liniments','Phosphate salts','Potassium Salt'

It overrides filters to allow: 'fludarabine', 'nilotinib', 'diamorphine', 'cocaine', 'lopinavir', 'flucloxacillin', 'isoniazid', 'hydrochlorothiazide', 'tolbutamide', 'streptomycin', 'ampicillin', 'phenobarbital', 'azithromycin', 'tenofovir', 'peramivir', 'artemisinin', 'lapatinib', 'rociletinib'

It specifically removes the following chemicals: 'Cholesterol','Carbon dioxide','Adenine','Guanine','Thymine','Uracil','Cytosine','Adenosine','Guanosine','5-Methyluridine','Uridine','Cytidine','Deoxyadenosine','Deoxyguanosine','Thymidine','Deoxyuridine','Deoxycytidine','Heparin','Hydrocortisone','Estradiol','Tretinoin','Testosterone','Progesterone','Melatonin'

It also removes chemicals without any product names or FDA labels. To supplement DrugBank and PharmGKB MeSH mappings, it tries exact matching against MeSH terms with the name of the chemical, and the name of the chemical + 'hydrochloride'.

## Identifying Sentences of Interest

The [findPGxSentences.py](https://github.com/jakelever/pgxmine/blob/master/findPGxSentences.py) finds star alleles using a regular expression and then parses the document to find sentences that mention at least one chemical and at least one variant (including protein/DNA modifictions, rsIDs, and star alleles). It filters out sentences by requiring that one of the below substrings is found in the sentence or that an rsID is found. These strings are from pgx_filter_terms.txt.

- associat
- response
- study
- studies
- significan
- survival
- compare
- observe
- interaction
- efficacy
- toxicity
- metabol
- resistan
- sensitiv
- marker
- correlate
- outcome
- susceptibility
- pharmaco
- predict
- efflux

## Mapping rsIDs to gene names using dbSNP

In order to show gene names, rsIDs were mapped to them using dbSNP. The [linkRSIDToGeneName.py](https://github.com/jakelever/pgxmine/blob/master/linkRSIDToGeneName.py) script parses the [VCF download](ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF) of dbSNP. It uses the PubTator Central annotations to check with rsIDs need to be saved so that the entirety of dbSNP is not saved.

## Mapping star alleles to rsIDs (where appropriate)

Some star alleles contain only a single mutation and are represented by an rsID (e.g. POR\*28 to rs1057868). Most star alleles are more complex haplotypes of a series of mutations in a genes. No mapping exists for all star alleles, though some exist for some gene familes (e.g. PharmVar). The [linkStarToRSID.py](https://github.com/jakelever/pgxmine/blob/master/linkStarToRSID.py) script does a regular expression search of the biomedical literature to find occurrences of star alleles right next to an rsID, e.g. "We study POR\*28 (rs1057868)" with the rsID in brackets. We then manually filtered this list to produce starRSMappings.tsv.

