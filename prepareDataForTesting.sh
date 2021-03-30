#!/bin/bash
#
#SBATCH --job-name=pubtator_test
#
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24576
#SBATCH --qos long

set -ex

echo "WARNING: Using miniature/empty drug and SNP metadata files for testing purposes. Don't use for main run!"

mkdir data

touch data/dbsnp_selected.tsv
echo '[
{
	"PharmGKB": "PA450085",
	"name": "Irinotecan",
	"DrugBank": "DB00762",
	"MeSH": [
	  "C051890",
	  "D000077146"
	],
	"isCancerDrug": true
},
{
	"PharmGKB": "PA451363",
	"name": "Simvastatin",
	"DrugBank": "DB00641",
	"MeSH": [
	  "D019821"
	],
	"isCancerDrug": false
}
]
' > data/selected_chemicals.json

echo -e "54658\tUGT1A1\tprotein-coding" > data/gene_names.tsv
echo -e "1576\tCYP3A4\tprotein-coding" >> data/gene_names.tsv
echo -e "1577\tCYP3A5\tprotein-coding" >> data/gene_names.tsv

# Unzip the annotated training data of pharmacogenomics relations
gunzip -c annotations.variant_other.bioc.xml.gz > data/annotations.variant_other.bioc.xml
gunzip -c annotations.variant_star_rs.bioc.xml.gz > data/annotations.variant_star_rs.bioc.xml

