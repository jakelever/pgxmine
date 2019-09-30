source('dependencies.R')

pharmGKBFilenames <- c('var_drug_ann.tsv','var_fa_ann.tsv','var_pheno_ann.tsv')
pharmGKB <- as.data.table(matrix(nrow=0,ncol=2))
pharmGKB_pmids <- c()
colnames(pharmGKB) <- c('Variant','Chemical')
pharmGKB_modifiedDate <- "0000-00-00"
for (pharmGKBFilename in pharmGKBFilenames) {
  pharmGKBFilename <- normalizePath(pharmGKBFilename)
  fileInfo <- file.info(pharmGKBFilename)
  tmpPharmGKB_modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]
  if (tmpPharmGKB_modifiedDate > pharmGKB_modifiedDate) {
    pharmGKB_modifiedDate <- tmpPharmGKB_modifiedDate
    paper.pharmgkbDate <- format(as.Date(fileInfo$mtime), format="%d %B %Y")
  }
  
  tempPharmGKB <- fread(pharmGKBFilename,sep='\t',header=T,stringsAsFactors=T,quote='')
  pharmGKB_pmids <- c(pharmGKB_pmids,tempPharmGKB$PMID)
  tempPharmGKB <- tempPharmGKB[,c('Variant','Chemical')]
  
  pharmGKB <- rbind(pharmGKB,tempPharmGKB)
}
pharmGKB_pmids <- unique(sort(pharmGKB_pmids))

#pharmGKB <- pharmGKB[grep(",",pharmGKB$Chemical,fixed=TRUE),]
pharmGKB$Chemical <- gsub('"','',pharmGKB$Chemical)


# Unroll PharmGKB
s <- strsplit(pharmGKB$Chemical, split = ",")
pharmGKB <- data.frame(Variant = rep(pharmGKB$Variant, sapply(s, length)), Chemical = unlist(s))
pharmGKB$Chemical <- trimws(pharmGKB$Chemical)

s <- strsplit(as.character(pharmGKB$Variant), split = ",")
pharmGKB <- data.frame(Chemical = rep(pharmGKB$Chemical, sapply(s, length)), Variant = unlist(s))
pharmGKB$Variant <- trimws(pharmGKB$Variant)

pharmGKB$Chemical_Name <- gsub("\\s*\\(\\S*\\)",'',pharmGKB$Chemical, perl=T)
pharmGKB$Chemical_ID <- gsub(")","",gsub(".*\\(",'',pharmGKB$Chemical, perl=T))
pharmGKB$Chemical_ID[grep("PA[0-9]*",pharmGKB$Chemical,perl=T,invert=T)] <- ''

pharmGKB2 <- pharmGKB
pharmGKB2$Variant <- gsub(":\\d\\d$","",pharmGKB2$Variant, perl=T)

pharmGKB3 <- pharmGKB
pharmGKB3$Variant <- gsub(":\\d\\d:\\d\\d$","",pharmGKB3$Variant, perl=T)

pharmGKB <- rbind(pharmGKB,pharmGKB2,pharmGKB3)
pharmGKB <- pharmGKB[!duplicated(pharmGKB),]
pharmGKB3$Variant <- gsub(":\\d\\d:\\d\\d$","",pharmGKB3$Variant, perl=T)

pharmGKB2 <- NA
pharmGKB3 <- NA

pharmGKB$Chemical_ID_Variant <- paste(pharmGKB$Chemical_ID,pharmGKB$Variant)
pharmGKB$Chemical_ID_Variant[pharmGKB$Chemical_ID==''] <- ''

sentencesFilename <- 'pgxmine_sentences.tsv'
collatedFilename <- 'pgxmine_collated.tsv'
collated <- fread(collatedFilename,sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')

sentences <- fread(sentencesFilename,sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')
collated$Chemical_ID_Variant <- paste(collated$chemical_pharmgkb_id,collated$variant_id)

# Strip out suballeles from star allelles (e.g. CYP3A5*3A -> CYP3A5*3)
starAlleles <- grep("*",pharmGKB$Variant,fixed=T)
pharmGKB$Chemical_ID_Variant[starAlleles] <- gsub("[A-Z]$","",pharmGKB$Chemical_ID_Variant[starAlleles])
starAlleles <- grep("*",collated$variant_id,fixed=T)
collated$Chemical_ID_Variant[starAlleles] <- gsub("[A-Z]$","",collated$Chemical_ID_Variant[starAlleles])

pharmGKBAssociations <- unique(pharmGKB$Chemical_ID_Variant)
pgxMineAssociations <- unique(collated$Chemical_ID_Variant)
pharmGKBAssociations <- pharmGKBAssociations[pharmGKBAssociations!='']
pgxMineAssociations <- pgxMineAssociations[pharmGKBAssociations!='']
pharmGKBAssociations <- pharmGKBAssociations[!is.na(pharmGKBAssociations)]
pgxMineAssociations <- pgxMineAssociations[!is.na(pgxMineAssociations)]

associatonComparison <- venn.diagram(
  x = list(PharmGKB=pharmGKBAssociations , PGxMine=pgxMineAssociations ),
  scaled=T,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
associationComparisonPlot <- gTree(children=associatonComparison)
associationComparisonPlot <- grid.arrange(associationComparisonPlot,bottom='Chemical/Variant Associations')

paper.pharmGKBPercentageOverlap <- round(100*sum(pgxMineAssociations %in% pharmGKBAssociations) / length(pgxMineAssociations),1)
paper.associationsNotInPharmGKB <- prettyNum(length(pgxMineAssociations) - sum(pgxMineAssociations %in% pharmGKBAssociations),big.mark=',')

pgxmine_PMIDs <- unique(sentences$pmid)

pmidComparisonPlot <- venn.diagram(
  x = list(PharmGKB=pharmGKB_pmids , PGxMine=pgxmine_PMIDs ),
  scaled=T,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
pmidComparisonPlot <- gTree(children=pmidComparisonPlot)
pmidComparisonPlot <- grid.arrange(pmidComparisonPlot,bottom='PubMed IDs')


fig_comparison <- arrangeGrob(associationComparisonPlot,pmidComparisonPlot,ncol=2)
grid.arrange(fig_comparison)


pharmGKB$Chemical_ID

notInPGMine_id <- pharmGKB$Chemical_ID[!(pharmGKB$Chemical_ID %in% collated$chemical_pharmgkb_id)]
notInPGMine <- sort(unique(pharmGKB[pharmGKB$Chemical_ID %in% notInPGMine_id,'Chemical_Name']))

#pmids_PMCAMC <- as.integer(scan("pmids/PMCAMC.txt", character(), quote = ""))
#pmids_PMCOA <- as.integer(scan("pmids/PMCOA.txt", character(), quote = ""))
#pmids_fulltext <- unique(c(pmids_PMCAMC,pmids_PMCOA))
#paper.pharmgkbPaperHasFullText <- prettyNum(sum(pharmGKB_pmids %in% pmids_fulltext),big.mark=',')
#paper.pharmgkbPaperNum <- prettyNum(length(pharmGKB_pmids),big.mark=',')
#paper.pharmgkbPaperHasFullTextPerc <- round(100*sum(pharmGKB_pmids %in% pmids_fulltext)/length(pharmGKB_pmids),1)

# Hard-coded as PMID data files are too big for GitHub
paper.pharmgkbPaperHasFullText <- "796"
paper.pharmgkbPaperNum <- "6,489"
paper.pharmgkbPaperHasFullTextPerc <- "12.3"
