#!/bin/bash
set -ex

rm -f example/aligned.bioc.xml example/sentences.bioc.xml example/kb.tsv example/mini_unfiltered.tsv example/mini_collated.tsv example/mini_sentences.tsv

# Align the PubTator Central extracted entities against the text sources to get offset positions for chemicals, variants
python align.py --inBioc example/input.bioc.xml --annotations <(zcat data/bioconcepts2pubtatorcentral.gz) --outBioc example/aligned.bioc.xml

# Parse and find sentences that mention a chemical, variant and likely a pharmacogenomic assocation (using filter terms)
python findPGxSentences.py --inBioc example/aligned.bioc.xml --filterTermsFile pgx_filter_terms.txt --outBioc example/sentences.bioc.xml

# Train relation classifiers (using the annotations* files as training data), filter for specific chemicals and apply the classifiers to extract associations and output with normalized genes, variants and chemicals
python createKB.py --trainingFiles annotations.variant_star_rs.bioc.xml,annotations.variant_other.bioc.xml --inBioC example/sentences.bioc.xml --selectedChemicals selected_chemicals.json --dbsnp dbsnp_selected.tsv --variantStopwords stopword_variants.txt --genes gene_names.tsv  --outKB example/kb.tsv

# Collate the output of createKB (which in a full run would be ~1800 files) and filter using the relation probability and collated by counting number of papers
python filterAndCollate.py --inData example --outUnfiltered example/mini_unfiltered.tsv --outCollated example/mini_collated.tsv --outSentences example/mini_sentences.tsv

