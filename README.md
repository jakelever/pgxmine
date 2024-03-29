# PGxMine

<p>
<a href="https://pgxmine.pharmgkb.org/">
   <img src="https://img.shields.io/badge/data-viewer-9e42f4.svg" />
</a>
<a href="https://doi.org/10.5281/zenodo.3360930">
   <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3360931.svg" />
</a>
</p>

This is the codebase for the PGxMine project that uses text-mining to identify papers for curation into [PharmGKB](https://www.pharmgkb.org). It is a Python3 project that makes use of the [Kindred](https://github.com/jakelever/kindred) relation classifier along with the [BioText](https://github.com/jakelever/biotext) project to manage the download of PubMed/PMC and alignment with [PubTator](https://www.ncbi.nlm.nih.gov/research/pubtator/).

# Viewing the Data

The data can be viewed through the [Shiny app](https://pgxmine.pharmgkb.org/). It can be downloaded as TSV files at [Zenodo](https://doi.org/10.5281/zenodo.3360930).

To run a local instance of the PGxmine viewer, the R Shiny code can be found in [shiny/](https://github.com/jakelever/pgxmine/tree/master/shiny) and installation instructions are found there too.

# Software Dependencies

This project depends on [Kindred](https://github.com/jakelever/kindred), [scispacy](https://allenai.github.io/scispacy/) and [snakemake](https://snakemake.github.io/). They can be installed by:

```
pip install -r requirements.txt
```

# Data Dependencies

This project uses a variety of data sources. 


A few need to be downloaded as below, apart from DrugBank which needs to be download manually.

- [MeSH](https://www.nlm.nih.gov/databases/download/mesh.html) (only needed to update the drug list)
- [DrugBank](https://www.drugbank.ca/releases/latest) (download manually as account is required and name it drugbank.xml)
- [PharmGKB](https://www.pharmgkb.org/downloads) (used for constructing the drug list and comparisons)

The [prepareData.sh](https://github.com/jakelever/pgxmine/blob/master/prepareData.sh) script downloads some of the data dependencies and runs some preprocessing to extract necessary data (such as gene name mappings). The commands that it runs are detailed below.

```
year=`date +"%Y"`

# Download MeSH, dbSNP, Entrez Gene metadata and pharmGKB drug info
sh downloadDataDependencies.sh

# Extract the gene names associated with rsIDs from dbSNP
python linkRSIDToGeneName.py --dbsnp <(zcat data/GCF_000001405.*.gz) --pubtator <(zcat data/bioconcepts2pubtatorcentral.gz) --outFile data/dbsnp_selected.tsv

# Create the drug list with mappings from MeSH IDs to PharmGKB IDs (with some filtering using DrugBank categories)
python createDrugList.py --meshC data/c$year.bin --meshD data/d$year.bin --drugbank drugbank.xml --pharmgkb data/drugs.tsv --outFile data/selected_chemicals.json

# Extract a mapping from Entrez Gene ID to name
zgrep -P "^9606\t" data/gene_info.gz | cut -f 2,3,10 -d $'\t' > data/gene_names.tsv

# Unzip the annotated training data of pharmacogenomics relations
gunzip -c annotations.variant_other.bioc.xml.gz > data/annotations.variant_other.bioc.xml
gunzip -c annotations.variant_star_rs.bioc.xml.gz > data/annotations.variant_star_rs.bioc.xml
```

# Example Run

There is an example input file in the [test\_data directory](https://github.com/jakelever/pgxmine/tree/master/test_data) which contains an PubMed abstract in BioC format. The [run\_example.sh script](https://github.com/jakelever/pgxmine/blob/master/run_example.sh) does a full run extracting chemical/variant associations and is shown below with comments. The final output is three files: mini\_unfiltered.tsv, mini\_collated.tsv, mini\_sentences.tsv. This is equivalent to the test run with snakemake shown below.

# Running with Snakemake

To run a small example of the pipeline using snakemake, run the command below. This runs on the data in the [test\_data directory](https://github.com/jakelever/pgxmine/tree/master/test_data). It is equivalent to the commands in the [run\_example.sh script](https://github.com/jakelever/pgxmine/blob/master/run_example.sh) which provides some comments on what each step does. Snakemake is useful for running on the larger datasets with the full run commands further down.

```
MODE=test snakemake --cores 1
```

To do a full run, you need set up a local instance of [BioText](https://github.com/jakelever/biotext) with the biocxml format. The command below will run Snakemake on the biotext. You must change BIOTEXT to point towards the biocxml directory in your local instance of BioText. The run will take a while and a cluster is recommended using snakemake's [cluster support](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

```
MODE=full BIOTEXT=/path/to/biotext/biocxml snakemake --cores 1
```

# Script Overview

Here is a summary of the main script files. The Snakefile manages the execution of these in the correct ordering.

## Main scripts

- **[findPGxSentences.py](https://github.com/jakelever/pgxmine/blob/master/findPGxSentences.py)**: Identify star alleles then find sentences that mention a chemical and variant
- **[getRelevantMeSH.py](https://github.com/jakelever/pgxmine/blob/master/getRelevantMeSH.py)**: Extracts MeSH terms related to age groups that is used by additional analysis
- **[createKB.py](https://github.com/jakelever/pgxmine/blob/master/createKB.py)**: Train and apply a relation classifier to extract pharmacogenomic chemical/variant associations
- **[filterAndCollate.py](https://github.com/jakelever/pgxmine/blob/master/filterAndCollate.py)**: Filter the results to reduce false positives and collate the associations
- **[utils/__init__.py](https://github.com/jakelever/pgxmine/blob/master/utils/__init__.py)**: Big functions for variant normalization and outputting the formatted sentences

## Other scripts

- **[createDrugList.py](https://github.com/jakelever/pgxmine/blob/master/createDrugList.py)**: Creates the list of drugs and drug mappings from MeSH IDs to PharmGKB IDs with some filtering by categories
- **[linkRSIDToGeneName.py](https://github.com/jakelever/pgxmine/blob/master/linkRSIDToGeneName.py)**: Extracts gene names from dbSNP associated with rsIDs
- **[linkStarToRSID.py](https://github.com/jakelever/pgxmine/blob/master/linkStarToRSID.py)**: Some rudimentary text mining to link star alleles with a specific rsID
- **[prepareForAnnotation.py](https://github.com/jakelever/pgxmine/blob/master/prepareForAnnotation.py)**: Select sentences and output to the standoff format to be annotated
- **[prCurve.py](https://github.com/jakelever/pgxmine/blob/master/prCurve.py)**: Calculate PR curves for the classifiers

# Paper

The paper can be recompiled using the dataset using Bookdown. All text and code for stats/figures are in the [paper/](https://github.com/jakelever/pgxmine/tree/master/paper) directory.

# Supplementary Materials

Supplementary materials for the manuscript are found in [supplementaryMaterials/](https://github.com/jakelever/pgxmine/blob/master/supplementaryMaterials/).

