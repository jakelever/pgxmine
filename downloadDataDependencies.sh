#!/bin/bash
set -ex

rm -fr data
mkdir -p data
cd data

year=`date +"%Y"`

wget ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/bioconcepts2pubtatorcentral.gz

wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz

wget ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/c$year.bin
wget ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/d$year.bin

wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz

wget https://api.pharmgkb.org/v1/download/file/data/drugs.zip
unzip drugs.zip

