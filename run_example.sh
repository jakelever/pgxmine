#!/bin/bash
set -ex

# Parse and find sentences that mention a chemical, variant and likely a pharmacogenomic assocation (using filter terms). Note that the input file already has chemicals and variants extracted from PubTator
python findPGxSentences.py --inBioc example/pubmed_26736037.bioc.xml --filterTermsFile pgx_filter_terms.txt --outBioc example/pubmed_26736037.sentences.bioc.xml

# Scan the documents for MeSH terms related to age which are provided as additional metadata later
python getRelevantMeSH.py --inBioc example/pubmed_26736037.bioc.xml --outJSONGZ example/pubmed_26736037.mesh.json.gz

# Train relation classifiers (using the annotations* files as training data), filter for specific chemicals and apply the classifiers to extract associations and output with normalized genes, variants and chemicals
python createKB.py --trainingFiles data/annotations.variant_star_rs.bioc.xml,data/annotations.variant_other.bioc.xml --inBioC example/pubmed_26736037.sentences.bioc.xml --selectedChemicals data/selected_chemicals.json --dbsnp data/dbsnp_selected.tsv --variantStopwords stopword_variants.txt --genes data/gene_names.tsv --relevantMeSH example/pubmed_26736037.mesh.json.gz  --outKB example/pubmed_26736037.kb.tsv

# Collate the output of createKB (which in a full run would be ~1800 files) and filter using the relation probability and collated by counting number of papers
python filterAndCollate.py --inData example --outUnfiltered example/mini_unfiltered.tsv --outCollated example/mini_collated.tsv --outSentences example/mini_sentences.tsv

