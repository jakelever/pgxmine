This describes the output files for the [PGxMine](https://github.com/jakelever/pgxmine) project. The code for this viewer is available in the PGxMine Github repo if you want to run it independently. Each file is a tab-delimited file with a header, no comments and no quoting.

You likely want **pgxmine\_collated.tsv** if you just want the list of chemical/variant assocations. If you want the supporting sentences, look at **pgxmine\_sentences.tsv**. You can use the *matching\_id* column to connect the two files. If you want to dig further and are okay with a higher false positive rate, look at **pgxmine\_unfiltered.tsv**.

**pgxmine\_collated.tsv:** This contains the chemical/variant associations with citation counts supporting them. It contains the normalized chemical, variant and where appropriate gene names with identifiers for PharmGKB, dbSNP and Entrez.

**pgxmine\_sentences.tsv:** This contains the supporting sentences for the chemical/variant associations in the collated file. Each row is a single supporting sentence for one association. This file contains information on the source publication (e.g. journal, publication date, etc), the actual sentence and the chemical/variant association extracted.

**pgxmine\_unfiltered.tsv:** This is the combined raw output of the createKB.py script across all of PubMed, Pubmed Central Open Access and PubMed Central Author Manuscript Collection. It contains every predicted relation with a prediction score above 0.5. So this may contain many false positives. Each row contain information on the publication (e.g. journal, publication date, etc) along with the sentence and the specific chemical/variant association.

