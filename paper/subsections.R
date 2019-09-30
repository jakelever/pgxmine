source('dependencies.R')

sentencesFilename <- 'pgxmine_sentences.tsv'

sentences <- fread(sentencesFilename,sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')


sentences$subsection <- as.character(sentences$subsection)
sentences[sentences$subsection=='conclusions','subsection'] <- 'conclusion'
sentences[sentences$subsection=='materials and methods','subsection'] <- 'methods'
sentences[sentences$subsection=='background','subsection'] <- 'introduction'
sentences[sentences$subsection=='materials and methods','subsection'] <- 'methods'
sentences[sentences$subsection=='study design','subsection'] <- 'methods'
sentences[sentences$subsection=='ethics statement','subsection'] <- 'methods'
sentences[sentences$subsection=='statistical analysis','subsection'] <- 'methods'
sentences[sentences$subsection=='data analysis','subsection'] <- 'methods'
sentences[sentences$subsection=='summary','subsection'] <- 'introduction'
sentences[sentences$subsection=='material and methods','subsection'] <- 'methods'
sentences[sentences$subsection=='method','subsection'] <- 'methods'
sentences[sentences$subsection=='statistical analyses','subsection'] <- 'methods'
sentences[sentences$subsection=='additional information','subsection'] <- 'methods'
sentences[sentences$subsection=='materials','subsection'] <- 'methods'
sentences[sentences$subsection=='data collection','subsection'] <- 'methods'
sentences[sentences$subsection=='case report','subsection'] <- 'methods'
sentences[sentences$subsection=='patients and method','subsection'] <- 'methods'
sentences[sentences$subsection=='patients and methods','subsection'] <- 'methods'
sentences[sentences$subsection=='limitations','subsection'] <- 'methods'
sentences[sentences$subsection=='results and discussion results','subsection'] <- 'results'
sentences[sentences$subsection=='results and discussion','subsection'] <- 'results'
sentences[sentences$subsection=='supplementary material','subsection'] <- 'methods'
sentences[sentences$subsection=='None','subsection'] <- 'Unable to identify'

sentences[sentences$subsection=='analysis','subsection'] <- 'methods'

sentences[sentences$subsection=='conflict of interest','subsection'] <- 'other'
sentences[sentences$subsection=='author contributions','subsection'] <- 'other'
sentences[sentences$subsection=='supporting information','subsection'] <- 'other'
sentences[sentences$subsection=='conflicts of interest','subsection'] <- 'other'
sentences[sentences$subsection=='abbreviations','subsection'] <- 'other'

sentences[sentences$subsection=='measures','subsection'] <- 'methods'
sentences[sentences$subsection=='statistics','subsection'] <- 'methods'
sentences[sentences$subsection=='statistical methods','subsection'] <- 'methods'
sentences[sentences$subsection=='participants','subsection'] <- 'methods'
#sentences[sentences$subsection!='abbreviations',]

sentences$subsection <- factor(as.character(sentences$subsection),levels=c('introduction','methods','results','discussion','conclusion','other','unknown','Unable to identify'))

subsectionCounts <- plyr::count(sentences[sentences$section=='article',c('subsection'),drop=F])
subsectionPlot <- barchart(freq ~ subsection, 
                           subsectionCounts,
                           ylab="Mentions",
                           ylim=c(0,1.1*max(subsectionCounts$freq)),
                           scales=list(x=list(rot=45)),
                           col="black")
grid.arrange(subsectionPlot)
