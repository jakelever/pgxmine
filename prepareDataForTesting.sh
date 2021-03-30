#!/bin/bash
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

echo "54658	UGT1A1	protein-coding" > data/gene_names.tsv
echo "1576	CYP3A4	protein-coding" >> data/gene_names.tsv
echo "1577	CYP3A5	protein-coding" >> data/gene_names.tsv

# Unzip the annotated training data of pharmacogenomics relations
gunzip -c annotations.variant_other.bioc.xml.gz > data/annotations.variant_other.bioc.xml
gunzip -c annotations.variant_star_rs.bioc.xml.gz > data/annotations.variant_star_rs.bioc.xml

