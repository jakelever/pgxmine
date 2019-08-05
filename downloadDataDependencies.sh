#!/bin/bash
set -ex

rm -fr data
mkdir -p data
cd data

wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.bgz

wget ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/c2019.bin
wget ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/d2019.bin

wget ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/bioconcepts2pubtatorcentral.gz

wget https://api.pharmgkb.org/v1/download/file/data/drugs.zip
unzip drugs.zip

