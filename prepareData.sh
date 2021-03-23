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

dummyDrugbank="<emptydrugbank></emptydrugbank>"

if [[ "$@" == "--useDummyDrugbank" ]]; then
	# For testing without DrugBank. Do not use for full run
	dummyDrugbankMD5=`echo "$dummyDrugbank" | md5sum`
	existingDrugbankMD5=`cat drugbank.xml | md5sum`
	if [ -f drugbank.xml && [ "$dummyDrugbankMD5" != "$existingDrugbankMD5" ] ]; then
		echo "ERROR: Drugbank.xml exists and contains non-dummy data. Cannot overwrite it with --useDummyDrugbank"
		exit 1
	fi
	echo "$dummyDrugbank" > drugbank.xml
elif [ ! -f drugbank.xml ]; then
	echo "ERROR: Couldn't find file drugbank.xml. Please download the latest version from https://www.drugbank.ca/releases and rename to drugbank.xml"
	exit 1
fi

dummyDrugbankMD5=`echo "$dummyDrugbank" | md5sum`
existingDrugbankMD5=`cat drugbank.xml | md5sum`
if [[ "$dummyDrugbankMD5" == "$existingDrugbankMD5" ]]; then
	echo "WARNING: Using a dummy version of DrugBank. You must download the actual DrugBank and save to drugbank.xml for full functionality"
fi

set -x
exit 0

# Download MeSH, dbSNP, Entrez Gene metadata and pharmGKB drug info
sh downloadDataDependencies.sh

# Extract the gene names associated with rsIDs from dbSNP
python linkRSIDToGeneName.py --dbsnp <(zcat data/GCF_000001405.25.gz) --pubtator <(zcat data/bioconcepts2pubtatorcentral.gz) --outFile data/dbsnp_selected.tsv

# Create the drug list with mappings from MeSH IDs to PharmGKB IDs (with some filtering using DrugBank categories)
python createDrugList.py --meshC data/c2019.bin --meshD data/d2019.bin --drugbank drugbank.xml --pharmgkb data/drugs.tsv --outFile data/selected_chemicals.json

# Extract a mapping from Entrez Gene ID to name
zgrep -P "^9606\t" data/gene_info.gz | cut -f 2,3,10 -d $'\t' > data/gene_names.tsv

# Unzip the annotated training data of pharmacogenomics relations
gunzip -c annotations.variant_other.bioc.xml.gz > data/annotations.variant_other.bioc.xml
gunzip -c annotations.variant_star_rs.bioc.xml.gz > data/annotations.variant_star_rs.bioc.xml

