#!/bin/bash
#
#SBATCH --job-name=pgxmine_data
#
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH -p rbaltman

set -e

if [ ! -f drugbank.xml ]; then
	echo "ERROR: Couldn't find file drugbank.xml. Please download the latest version from https://www.drugbank.ca/releases and rename to drugbank.xml"
	exit 1
fi

set -x

year=`date +"%Y"`

# Download MeSH, dbSNP, Entrez Gene metadata and pharmGKB drug info
sh downloadDataDependencies.sh

# Extract the gene names associated with rsIDs from dbSNP
python linkRSIDToGeneName.py --dbsnp <(zcat data/GCF_000001405.25.gz) --pubtator <(zcat data/bioconcepts2pubtatorcentral.gz) --outFile data/dbsnp_selected.tsv

# Create the drug list with mappings from MeSH IDs to PharmGKB IDs (with some filtering using DrugBank categories)
python createDrugList.py --meshC data/c$year.bin --meshD data/d$year.bin --drugbank drugbank.xml --pharmgkb data/drugs.tsv --outFile data/selected_chemicals.json

# Extract a mapping from Entrez Gene ID to name
zgrep -P "^9606\t" data/gene_info.gz | cut -f 2,3,10 -d $'\t' > data/gene_names.tsv

# Unzip the annotated training data of pharmacogenomics relations
gunzip -c annotations.variant_other.bioc.xml.gz > data/annotations.variant_other.bioc.xml
gunzip -c annotations.variant_star_rs.bioc.xml.gz > data/annotations.variant_star_rs.bioc.xml

