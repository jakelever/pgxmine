#!/bin/bash
#
#SBATCH --job-name=pubtator_test
#
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24576
#SBATCH --qos long

set -e

if [ ! -f drugbank.xml ]; then
	echo "ERROR: Couldn't find file drugbank.xml. Please download the latest version from https://www.drugbank.ca/releases and rename to drugbank.xml"
	exit 1
fi

set -x

# Download PubTator Central, MeSH, dbSNP, Entrez Gene metadata and pharmGKB drug info
sh downloadDataDependencies.sh

# Extract the gene names associated with rsIDs from dbSNP
python linkRSIDToGeneName.py --dbsnp <(zcat data/GCF_000001405.25.gz) --pubtator <(zcat data/bioconcepts2pubtatorcentral.gz) --outFile dbsnp_selected.tsv

# Create the drug list with mappings from MeSH IDs to PharmGKB IDs (with some filtering using DrugBank categories)
python createDrugList.py --meshC data/c2019.bin --meshD data/d2019.bin --drugbank drugbank.xml --pharmgkb data/drugs.tsv --outFile selected_chemicals.json

# Extract a mapping from Entrez Gene ID to name
zgrep -P "^9606\t" data/gene_info.gz | cut -f 2,3,10 -d $'\t' > gene_names.tsv

# Unzip the annotated training data of pharmacogenomics relations
gunzip -c annotations.variant_other.bioc.xml.gz > annotations.variant_other.bioc.xml
gunzip -c annotations.variant_star_rs.bioc.xml.gz > annotations.variant_star_rs.bioc.xml
