## Example document for mini run

This directory contains an example input file to the pipeline which is run with the [run\_example.sh](https://github.com/jakelever/pgxmine/blob/master/run_example.sh) script. The input file contains a [single abstract](https://pubmed.ncbi.nlm.nih.gov/26736037/) which has had the PubTator entities matched in its text. This functionality is now part of the [BioText](https://github.com/jakelever/biotext) project.

Specifically, to create this single test file, the commands below were used:
```
# Set $BIOTEXT to the directory with the BioText project code.
# Set $EMAIL to your email (as the NCBI API requires your email)
# You also need to download the PubTator annotations file which can be done with 'snakemake --cores 1 pubtator_downloaded.flag'

# This fetches the abstract from the PubMed database
python $BIOTEXT/src/convertEUtils.py --database pubmed --identifiers 26736037 --email $EMAIL --o pubmed_26736037.noentities.bioc.xml --oFormat biocxml

# This uses the PubTator dataset to identify mentions of biomedical entities
python $BIOTEXT/src/alignWithPubtator.py --inBioc pubmed_26736037.noentities.bioc.xml --annotations <(zcat $BIOTEXT/bioconcepts2pubtatorcentral.gz) --outBioc pubmed_26736037.bioc.xml
```

