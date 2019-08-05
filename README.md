# PGxMine

This is the codebase for the PGxMine project to using text-mining to identify papers for curation into PharmGKB. It is a Python3 project that makes use of the Kindred relation classifier along with the PubRunner project to run tools across PubMed and accessible PubMed Central.

# Software Dependencies

This project dependends on Kindred, scispacy and PubRunner. They can be installed as below:

```
pip install kindred pubrunner scispacy

pip install scispacy
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.2.0/en_core_sci_md-0.2.0.tar.gz
```

# Data Dependencies

This project uses a variety of data sources. The downloadPubtatorCentralAndMeSH.sh script will download two and PubRunner will manage the others, apart from DrugBank which needs to be download manually.

- PubMed and accessible PubMed Central (downloaded by PubRunner)
- PubTator Central
- MeSH (only needed to update the drug list)
- DrugBank (download manually as account is required and name it full_database.xml)
- PharmGKB (used for constructing the drug list and comparisons)

The following commands prepare the needed data (after you have downloaded the DrugBank XML file manually).

```
# Download PubTator Central, MeSH and dbSNP
sh downloadDataDependencies.sh

# Extract the gene names associated with rsIDs from dbSNP
python linkRSIDToGeneName.py --dbsnp <(zcat data/GCF_000001405.25.bgz) --pubtator <(zcat data/bioconcepts2pubtatorcentral.gz) --outFile gene_names.tsv

# Create the drug list with mappings from MeSH IDs to PharmGKB IDs (with some filtering using DrugBank categories)
python createDrugList.py --meshC data/c2019.bin --meshD data/d2019.bin --drugbank full_database.xml --pharmgkb data/drugs.tsv --outFile selected_chemicals.json
```

# Example Run

# Full Run

# Viewer

# Supplementary Materials

In the FILE, there are additional details of methods that didn't fit into the paper.

