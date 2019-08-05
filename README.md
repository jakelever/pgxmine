# PGxMine

This is the codebase for the PGxMine project to using text-mining to identify papers for curation into [PharmGKB](https://www.pharmgkb.org). It is a Python3 project that makes use of the [Kindred](https://github.com/jakelever/kindred) relation classifier along with the [PubRunner](https://github.com/jakelever/pubrunner) project to run tools across PubMed and accessible PubMed Central.

# Viewing the Data

The data can be viewed through the [Shiny app](http://bionlp.bcgsc.ca/pgxmine/).

To run a local instance of the PGxmine viewer, the R Shiny code can be found in [shiny/](https://github.com/jakelever/pgxmine/tree/master/shiny) and installation instructions are found there too.

# Software Dependencies

This project depends on [Kindred](https://github.com/jakelever/kindred), [scispacy](https://allenai.github.io/scispacy/) and [PubRunner](https://github.com/jakelever/pubrunner). They can be installed as below:

```
pip install kindred pubrunner scispacy

pip install scispacy
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.2.0/en_core_sci_md-0.2.0.tar.gz
```

# Data Dependencies

This project uses a variety of data sources. The [downloadDataDependencies.sh](https://github.com/jakelever/pgxmine/blob/master/downloadDataDependencies.sh) script will download two and PubRunner will manage the others, apart from DrugBank which needs to be download manually.

- [PubMed](https://www.nlm.nih.gov/databases/download/pubmed_medline.html) and accessible [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/tools/ftp/) (downloaded by PubRunner)
- [PubTator Central](https://www.ncbi.nlm.nih.gov/research/pubtator/)
- [MeSH](https://www.nlm.nih.gov/databases/download/mesh.html) (only needed to update the drug list)
- [DrugBank](https://www.drugbank.ca/releases/latest) (download manually as account is required and name it full\_database.xml)
- [PharmGKB](https://www.pharmgkb.org/downloads) (used for constructing the drug list and comparisons)

The following commands prepare the needed data (after you have downloaded the DrugBank XML file manually).

```
# Download PubTator Central, MeSH and dbSNP
sh downloadDataDependencies.sh

# Extract the gene names associated with rsIDs from dbSNP
python linkRSIDToGeneName.py --dbsnp <(zcat data/GCF_000001405.25.bgz) --pubtator <(zcat data/bioconcepts2pubtatorcentral.gz) --outFile gene_names.tsv

# Create the drug list with mappings from MeSH IDs to PharmGKB IDs (with some filtering using DrugBank categories)
python createDrugList.py --meshC data/c2019.bin --meshD data/d2019.bin --drugbank full_database.xml --pharmgkb data/drugs.tsv --outFile selected_chemicals.json
```

# Example Run

There is an example input file in the example directory which contains a couple PubMed abstracts in BioC format. The run\_example.sh script does a full run extracting chemical/variant associations and is shown below with comments. The final output is three files: mini\_unfiltered.tsv, mini\_collated.tsv, mini\_sentences.tsv.

```
# Align the PubTator Central extracted entities against the text sources to get offset positions for chemicals, variants
python align.py --inBioc example/input.bioc.xml --annotations <(zcat data/bioconcepts2pubtatorcentral.gz) --outBioc example/aligned.bioc.xml

# Parse and find sentences that mention a chemical, variant and likely a pharmacogenomic assocation (using filter terms)
python findPGxSentences.py --inBioc example/aligned.bioc.xml --filterTermsFile pgx_filter_terms.txt --outBioc example/sentences.bioc.xml

# Train relation classifiers (using the annotations* files as training data), filter for specific chemicals and apply the classifiers to extract associations and output with normalized genes, variants and chemicals
python createKB.py --trainingFiles annotations.variant_star_rs.bioc.xml,annotations.variant_other.bioc.xml --inBioC example/sentences.bioc.xml --selectedChemicals selected_chemicals.json --dbsnp dbsnp_selected.tsv --variantStopwords stopword_variants.txt --genes gene_names.tsv  --outKB example/kb.tsv

# Collate the output of createKB (which in a full run would be ~1800 files) and filter using the relation probability and collated by counting number of papers
python filterAndCollate.py --inData example --outUnfiltered example/mini_unfiltered.tsv --outCollated example/mini_collated.tsv --outSentences example/mini_sentences.tsv
```

# Full Run

To do a full run, you would need to use PubRunner. It would then manage the download and format conversion of PubMed, PubMed Central Open Access subset and PubMed Central Author Manuscript Collection. The command to do it is below

```
pubrunner .
```

This will take a long time. Setting up PubRunner with a cluster is recommended. A test run is below.

```
pubrunner --test .
```

# Script Overview

Here is a summary of the main script files. The pubrunner.yml file is the master script for PubRunner and lists the resources and script usage to actually run the project.

## Setup scripts

- **createDrugList.py**: Creates the list of drugs and drug mappings from MeSH IDs to PharmGKB IDs with some filtering by categories
- **linkRSIDToGeneName.py**: Extracts gene names from dbSNP associated with rsIDs
- **linkStarToRSID.py**: Some rudimentary text mining to link star alleles with a specific rsID
- **prCurve.py**: Calculate PR curves for the classifiers

## Main scripts

- **align.py**: Align PubTator Central entities against abstracts and full-text papers
- **findPGxSentences.py**: Identify star alleles then find sentences that mention a chemical and variant
- **createKB.py**: Train and apply a relation classifier to extract pharmacogenomic chemical/variant associations
- **filterAndCollate.py**: Filter the results to reduce false positives and collate the associations

# Paper

The paper can be recompiled using the dataset using Bookdown. All text and code for stats/figures are in the [paper/](https://github.com/jakelever/pgxmine/tree/master/paper) directory.

# Supplementary Materials

Supplementary materials for the manuscript are found in [supplementary\_materials.md](https://github.com/jakelever/pgxmine/blob/master/supplementary_materials.md).

