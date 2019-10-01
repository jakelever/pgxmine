source('dependencies.R')

sentences <- fread('pgxmine_sentences.tsv',sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')
collated <- fread('pgxmine_collated.tsv',sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')

paper.mentionCount <- prettyNum(nrow(sentences),big.mark=',')
paper.pmidCount <- prettyNum(length(unique(sentences$pmid)),big.mark=',')
paper.sentenceCount <- prettyNum(length(unique(sentences$sentence)),big.mark=',')
paper.associationCount <- prettyNum(nrow(collated),big.mark=',')

paper.chemicalCount <- prettyNum(length(unique(collated$chemical_normalized)),big.mark=',')
paper.variantCount <- prettyNum(length(unique(collated$variant_id)),big.mark=',')

paper.titleOrAbstract <- round(100*sum(sentences$section %in% c('title','abstract')) / nrow(sentences),1)

#genes <- gene_name = unlist(strsplit(as.character(collated$gene_names), split = ","))
#paper.geneCount <- prettyNum(length(unique(collated$gene_names[collated$gene_names!=''])))

paper.normalizedToRSID <- round(100*length(grep("^rs",sentences$variant_id)) / nrow(sentences),1)
paper.starAlleles <- round(100*length(grep("^\\*",sentences$variant_normalized)) / nrow(sentences),1)

paper.singlePapers <- round(100*sum(collated$citation_count==1) / nrow(collated),1)

topAssoc <- collated[order(collated$citation_count,decreasing=T)[1],]
paper.topAssoc_chemical <- topAssoc$chemical_normalized
paper.topAssoc_variant <- topAssoc$variant_id
paper.topAssoc_paperCount <- topAssoc$citation_count

variantCounts <- plyr::count(paste(sentences$gene_names, sentences$variant_normalized))
variantCounts <- variantCounts[order(variantCounts$freq,decreasing=T),]
variantCounts <- variantCounts[1:10,]
variantCounts$x <- gsub('EGFR,EGFR-AS1','EGFR',variantCounts$x)
variantCounts$x <- factor(variantCounts$x, levels=as.character(variantCounts$x))
variantPlot <- barchart(freq ~ x, 
         variantCounts, 
         col='black', 
         xlab='Variant', 
         ylab='# of Extracted Associations',
         scales=list(x=list(rot=45)))

chemicalCounts <- plyr::count(sentences$chemical_normalized)
chemicalCounts <- chemicalCounts[order(chemicalCounts$freq,decreasing=T),]
chemicalCounts <- chemicalCounts[1:10,]
chemicalCounts$x <- factor(chemicalCounts$x, levels=as.character(chemicalCounts$x))
chemicalPlot <- barchart(freq ~ x, 
         chemicalCounts, 
         col='black', 
         xlab='Chemical', 
         ylab='# of Extracted Associations',
         scales=list(x=list(rot=45)))

journalCounts <- plyr::count(tolower(sentences$journal_short))
journalCounts <- journalCounts[order(journalCounts$freq,decreasing=T),]
journalCounts <- journalCounts[1:10,]
journalCounts$x <- factor(journalCounts$x, levels=as.character(journalCounts$x))
journalPlot <- barchart(freq ~ x, 
                         journalCounts, 
                         col='black', 
                         xlab='Journal', 
                         ylab='# of Extracted Associations',
                         scales=list(x=list(rot=45)))

fig_variantChemicalJournals <- arrangeGrob(chemicalPlot,variantPlot,journalPlot,nrow=1)
grid.arrange(fig_variantChemicalJournals)

prCurve_other <- read.table('prcurve.variant_other.tsv',header=T)
prCurve_other <- prCurve_other[order(prCurve_other$precision),]
prCurve_other <- prCurve_other[order(prCurve_other$recall),]
prCurve_other$source <- 'DNA & Protein Modifications'

prCurve_star_rs <- read.table('prcurve.variant_star_rs.tsv',header=T)
prCurve_star_rs <- prCurve_star_rs[order(prCurve_star_rs$precision),]
prCurve_star_rs <- prCurve_star_rs[order(prCurve_star_rs$recall),]
prCurve_star_rs$source <- 'Star Alleles & RS IDs'

myColours <- c(brewer.pal(3,"Dark2"),"#000000")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black"),
  strip.background=list(col="lightgrey")
)

prCurve <- rbind(prCurve_other,prCurve_star_rs)
prCurve <- prCurve[order(prCurve$threshold),]
thresholdPlot <- xyplot(precision + recall ~ threshold | source, 
       prCurve,
       xlim=c(-0.1,1.1),
       ylim=c(-0.1,1.1),
       xlab='Threshold',
       ylab='Precision & Recall',
       lwd=2,
       par.settings = my.settings,
       auto.key=list(space="top", columns=2, 
                     points=FALSE, rectangles=TRUE),
       panel = function(x, y,  ...) {
         panel.abline(v = 0.75, col='black', lwd=1)
         panel.xyplot(x, y, ...)
       },
       type='l')

prCurve <- prCurve[order(prCurve$precision),]
prCurve <- prCurve[order(prCurve$recall),]
prcurvesPlot <- xyplot(precision ~ recall | source, 
       prCurve,
       xlim=c(-0.1,1.1),
       ylim=c(-0.1,1.1),
       xlab='Recall',
       ylab='Precision',
       par.settings = my.settings,
       lwd=2,
       col='black',
       type='l')

fig_prcurves <- arrangeGrob(prcurvesPlot,thresholdPlot,ncol=2)
grid.arrange(fig_prcurves)

atThreshold1 = prCurve[prCurve$source=='DNA & Protein Modifications' & prCurve$threshold==.75,]
atThreshold2 = prCurve[prCurve$source=='Star Alleles & RS IDs' & prCurve$threshold==.75,]
paper.precision1 <- round(100*atThreshold1$precision,1)
paper.precision2 <- round(100*atThreshold2$precision,1)
paper.recall1 <- round(100*atThreshold1$recall,1)
paper.recall2 <- round(100*atThreshold2$recall,1)

paper.avg_precision <- round((paper.precision1+paper.precision2)/2,1)
paper.avg_recall <- round((paper.recall1+paper.recall2)/2,1)

inTitle <- sum(sentences$section=='title')
inAbstract <- sum(sentences$section=='abstract')

paper.mentionsInFullText <- prettyNum(nrow(sentences) - inTitle - inAbstract,big.mark=',')
paper.mentionsInFullTextPerc <- round(100*(nrow(sentences) - inTitle - inAbstract)/nrow(sentences),1)
