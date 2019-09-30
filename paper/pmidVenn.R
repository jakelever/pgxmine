
source('dependencies.R')

pmids_PMCAMC <- as.integer(scan("pmids/PMCAMC.txt", character(), quote = ""))
pmids_PMCOA <- as.integer(scan("pmids/PMCOA.txt", character(), quote = ""))
pmids_PUBMED <- as.integer(scan("pmids/PUBMED.txt", character(), quote = ""))
pmids_PTC <- as.integer(scan("pmids/pubtator_central.txt", character(), quote = ""))

pmidCoveragePlot <- venn.diagram(
  x = list(PMCAMC=pmids_PMCAMC , PMCOA=pmids_PMCOA, PubMed=pmids_PUBMED,"PubTator Central"=pmids_PTC ),
  scaled=T,
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
pmidCoveragePlot <- gTree(children=pmidCoveragePlot)
grid.arrange(pmidCoveragePlot)

