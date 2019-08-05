#!/bin/bash
set -ex

rm -f example/aligned.bioc.xml example/sentences.bioc.xml example/kb.tsv example/mini_unfiltered.tsv example/mini_collated.tsv example/mini_sentences.tsv

python align.py --inBioc example/input.bioc.xml --annotations <(zcat data/bioconcepts2pubtatorcentral.gz) --outBioc example/aligned.bioc.xml
python findPGxSentences.py --inBioc example/aligned.bioc.xml --filterTermsFile pgx_filter_terms.txt --outBioc example/sentences.bioc.xml
python createKB.py --trainingFiles annotations.variant_star_rs.bioc.xml,annotations.variant_other.bioc.xml --inBioC example/sentences.bioc.xml --selectedChemicals selected_chemicals.json --dbsnp dbsnp_selected.tsv --variantStopwords stopword_variants.txt --genes gene_names.tsv  --outKB example/kb.tsv
python filterAndCollate.py --inData example --outUnfiltered example/mini_unfiltered.tsv --outCollated example/mini_collated.tsv --outSentences example/mini_sentences.tsv

