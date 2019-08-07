library(plotly)
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(stringr)

# Weird hack as R sometimes "forgets" its working directory
wd <- setwd(".")
setwd(wd)

sentencesFilename <- 'pgxmine_sentences.tsv'
collatedFilename <- 'pgxmine_collated.tsv'

sentencesFilename <- normalizePath(sentencesFilename)
collatedFilename <- normalizePath(collatedFilename)
fileInfo <- file.info(sentencesFilename)
pgxmine_modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]

sentences <- fread(sentencesFilename,sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')
collated <- fread(collatedFilename,sep='\t',header=T,stringsAsFactors=T,quote='', encoding = 'UTF-8')

sentences <- sentences[order(sentences$year,sentences$month,sentences$day,decreasing=T),]

pharmGKBFilenames <- c('var_drug_ann.tsv','var_fa_ann.tsv','var_pheno_ann.tsv')
pharmGKB <- as.data.table(matrix(nrow=0,ncol=2))
pharmGKB_pmids <- c()
colnames(pharmGKB) <- c('Variant','Chemical')
pharmGKB_modifiedDates <- c()
for (pharmGKBFilename in pharmGKBFilenames) {
  pharmGKBFilename <- normalizePath(pharmGKBFilename)
  fileInfo <- file.info(pharmGKBFilename)
  tmpPharmGKB_modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]
  pharmGKB_modifiedDates <- c(pharmGKB_modifiedDates,tmpPharmGKB_modifiedDate)
  
  tempPharmGKB <- fread(pharmGKBFilename,sep='\t',header=T,stringsAsFactors=T,quote='')
  pharmGKB_pmids <- c(pharmGKB_pmids,tempPharmGKB$PMID)
  tempPharmGKB <- tempPharmGKB[,c('Variant','Chemical')]
  
  pharmGKB <- rbind(pharmGKB,tempPharmGKB)
}
pharmGKB_pmids <- unique(sort(pharmGKB_pmids))

# Get the most recent modified date of the PharmGKB files loaded
pharmGKB_modifiedDate <- sort(pharmGKB_modifiedDates,decreasing=T)[1]

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

pharmGKB$Variant <- gsub("*0","*", pharmGKB$Variant, fixed=T)
collated$variant_id_edited <- gsub("*0","*", collated$variant_id, fixed=T)

#collated$variant_id_edited <- as.character(collated$variant_id_edited)
#pharmGKB$Variant <- as.character(pharmGKB$Variant)

starRSMapping <- read.table('starRSMappings.tsv',sep='\t',header=T)
starRSMapping$geneAndAllele <- paste(starRSMapping$gene,starRSMapping$allele,sep='*')
for(i in 1:nrow(starRSMapping)) {
  collated$variant_id_edited[collated$variant_id_edited==starRSMapping$geneAndAllele[i]] <- as.character(starRSMapping$rsid[i])
  pharmGKB$Variant[pharmGKB$Variant==starRSMapping$geneAndAllele[i]] <- as.character(starRSMapping$rsid[i])
}


pharmGKB$Chemical_ID_Variant <- paste(pharmGKB$Chemical_ID,pharmGKB$Variant)
pharmGKB$Chemical_ID_Variant[pharmGKB$Chemical_ID==''] <- ''

collated$Chemical_ID_Variant <- paste(collated$chemical_pharmgkb_id,collated$variant_id_edited)

# Strip out suballeles from star allelles (e.g. CYP3A5*3A -> CYP3A5*3)
starAlleles <- grep("*",pharmGKB$Variant,fixed=T)
pharmGKB$Chemical_ID_Variant[starAlleles] <- gsub("[A-Z]$","",pharmGKB$Chemical_ID_Variant[starAlleles])
starAlleles <- grep("*",collated$variant_id,fixed=T)
collated$Chemical_ID_Variant[starAlleles] <- gsub("[A-Z]$","",collated$Chemical_ID_Variant[starAlleles])


#collated$Chemical_ID_Variant <- gsub('*57:01:01','*5701',collated$Chemical_ID_Variant,fixed=T)
#collated$Chemical_ID_Variant <- gsub('*57:01','*5701',collated$Chemical_ID_Variant,fixed=T)


sentences$pmid_in_pharmgkb <- sentences$pmid %in% pharmGKB_pmids

#chemicalIgnoreList <- c('statin')
#sentences <- sentences[!(tolower(sentences$chemical_intext) %in% chemicalIgnoreList),]

collatedNoNA <- is.na(collated$gene_names)
s <- strsplit(as.character(collated$gene_names), split = ",")
genesToMatchingID <- data.frame(matching_id = rep(collated$matching_id, sapply(s, length)), gene_name = unlist(s))
genesToMatchingID$gene_name <- trimws(genesToMatchingID$gene_name)

chemicals <- sort(unique(as.character(collated$chemical_normalized)))
genes <- sort(unique(as.character(genesToMatchingID$gene_name)))
variants <- sort(unique(as.character(collated$variant_normalized)))

variantTypes <- sort(unique(as.character(collated$variant_type)))

collated$chemical_in_pharmgkb <- collated$chemical_pharmgkb_id %in% pharmGKB$Chemical_ID

collated$variant_in_pharmgkb <- collated$variant_id_edited %in% pharmGKB$Variant
collated$association_in_pharmgkb <- collated$Chemical_ID_Variant %in% pharmGKB$Chemical_ID_Variant

collated$chemical_in_pharmgkb[is.na(collated$chemical_pharmgkb_id)] <- NA
collated$association_in_pharmgkb[is.na(collated$chemical_pharmgkb_id)] <- NA
collated$variant_in_pharmgkb[collated$variant_id==''] <- NA
collated$association_in_pharmgkb[collated$variant_id==''] <- NA

associationsWithAnnotatedPaper <- unique(sentences$matching_id[sentences$pmid_in_pharmgkb])
collated$in_pharmgkb_paper <- F
collated[collated$matching_id %in% associationsWithAnnotatedPaper,'in_pharmgkb_paper'] <- T

sentences$pubmed_link <- paste("<a target=\"_blank\" href='https://www.ncbi.nlm.nih.gov/pubmed/", sentences$pmid, "'>", sentences$pmid, "</a>", sep='')

associationTableExplanation <- "<br /><b>Association Table:</b><br />The table below shows a list of extracted chemical-variant associations. Click on one and examine the table further below to see the individual sentences. Variants and chemicals are normalized back to PharmGKB entities. Blank values show where this has not been possible.<br /><br />"

paperSentenceTableExplanation <- "<br /><br /><br /><b>Sentence Table:</b><br />Select a chemical-variant association in the table above to see sentences and publication information<br /><br />"

lastModifiedText = paste("PGxMine updated on ",pgxmine_modifiedDate, " Comparing against PharmGKB last downloaded on ", pharmGKB_modifiedDate, sep="")



# Some cleanup
collated[,Chemical_ID_Variant:=NULL]


# https://stackoverflow.com/questions/27145572/dynamically-selecting-over-multiple-tabs-in-a-shiny-app
ui <- function(req) {
  fluidPage(
    tags$head(
      includeHTML("google-analytics.js"),
      tags$style(".rightAlign{float:right; margin-left:5px; margin-bottom: 20px;}")
    ),
    HTML('<script>
       function annotatePositive(id) {
           Shiny.onInputChange("annotate_positive_input", id);
           Shiny.onInputChange("trigger", Math.random());
           return false;
       };
       function annotateNegative(id) {
           Shiny.onInputChange("annotate_negative_input", id);  
           Shiny.onInputChange("trigger", Math.random());
           return false;
       };
         </script>'),         
    titlePanel("",windowTitle="PGxMine"),
    helpText(includeHTML("header.html")),
    #verbatimTextOutput("click"),
    #verbatimTextOutput("moo"),
    sidebarPanel(
      #sliderInput("threshold_input", "Classifier Threshold", min = 0.5, max = 1.0, value = 0.9, step=0.05),
      selectizeInput(inputId="chemical_input", 
                     label=p("Chemical",actionLink("chemical_clear", " (Clear)", style='font-size:70%')), 
                     choices=c('', chemicals), 
                     selected = '', 
                     multiple = FALSE, 
                     options = list(maxOptions = 2*length(chemicals))),
      selectizeInput(inputId = "gene_input", 
                     label=p("Gene",actionLink("gene_clear", " (Clear)", style='font-size:70%')), 
                     choices = c('', genes), 
                     selected = '', 
                     multiple = FALSE, 
                     options = list(maxOptions = 2*length(genes))),
      selectizeInput(inputId = "variant_input", 
                     label=p("Variant",actionLink("variant_clear", " (Clear)", style='font-size:70%')), 
                     choices = c('', variants), 
                     selected = '',
                     multiple = FALSE, 
                     options = list(maxOptions = 2*length(variants))),
      
      checkboxGroupInput("chemical_in_pharmgkb_input", "Chemical in PharmGKB", c('Yes', 'No')),
      checkboxGroupInput("polyinpharmgkb_input", "Variant in PharmGKB", c('Yes', 'No')),
      checkboxGroupInput("association_in_pharmgkb_input", "Association in PharmGKB", c('Yes', 'No')),
      checkboxGroupInput("in_pharmgkb_paper_input", "In PharmGKB Curated Paper", c('Yes', 'No')),
      
      checkboxGroupInput("varianttype_input", "Variant Type in Text", variantTypes),
      
      
      width=2
      #verbatimTextOutput("gene_text")
    ),
    mainPanel(
      HTML("<p></p>"),
      splitLayout(cellWidths = c("31%", "31%", "31%"), 
                  plotlyOutput("chemical_piechart"),
                  plotlyOutput("gene_piechart"),
                  plotlyOutput("variant_piechart")),
      
      HTML(associationTableExplanation),
      
      downloadButton("download_collated_all", label = "Download All", class='rightAlign'),
      downloadButton("download_collated_shown", label = "Download Shown", class='rightAlign'),
      
      DT::dataTableOutput("data_table"),
      
      HTML(paperSentenceTableExplanation),
      
      downloadButton("download_sentences_all", label = "Download All Sentences", class='rightAlign'),
      downloadButton("download_sentences_above", label = "Download Sentences for Chemical-Variants Above", class='rightAlign'),
      downloadButton("download_sentences_shown", label = "Download Sentences for Selected Chemical-Variant", class='rightAlign'),
      
      DT::dataTableOutput("sentences_table"),
      
      helpText(lastModifiedText)
    )
                  
                
    
    
  )
}

input <- data.frame(chemical_input='Oxytocin',gene_input='EGFR',variant_input='p.T790M')

server <- function(input, output, session) {
  
  sentenceData <- reactive({
    
  })
  
  tableData <- reactive({
    
    selected <- rep(TRUE,nrow(collated))
    
    if (!is.null(input$chemical_in_pharmgkb_input) && length(input$chemical_in_pharmgkb_input)==1) {
      if (input$chemical_in_pharmgkb_input == 'Yes') {
        selected <- selected & collated$chemical_in_pharmgkb==T
      } else if (input$chemical_in_pharmgkb_input == 'No') {
        selected <- selected & collated$chemical_in_pharmgkb==F
      }
    }
    
    
    if (!is.null(input$polyinpharmgkb_input)) {
      matching = c()
      if ('Yes' %in% input$polyinpharmgkb_input) {
        matching <- c(matching,T)
      }
      if ('No' %in% input$polyinpharmgkb_input) {
        matching <- c(matching,F)
      }
      selected <- selected & collated$variant_in_pharmgkb %in% matching
    }
    
    
    if (!is.null(input$association_in_pharmgkb_input)) {
      matching = c()
      if ('Yes' %in% input$association_in_pharmgkb_input) {
        matching <- c(matching,T)
      }
      if ('No' %in% input$association_in_pharmgkb_input) {
        matching <- c(matching,F)
      }
      selected <- selected & collated$association_in_pharmgkb %in% matching
    }
    
    
    if (!is.null(input$in_pharmgkb_paper_input)) {
      matching = c()
      if ('Yes' %in% input$in_pharmgkb_paper_input) {
        matching <- c(matching,T)
      }
      if ('No' %in% input$in_pharmgkb_paper_input) {
        matching <- c(matching,F)
      }
      selected <- selected & collated$in_pharmgkb_paper %in% matching
    }
    
    
    
    if (!is.null(input$varianttype_input) && length(input$varianttype_input)>0) {
      selected <- selected & collated$variant_type %in% input$varianttype_input
    }
    
    if (input$chemical_input!='') {
      selected <- selected & collated$chemical_normalized==input$chemical_input
    }
    
    if (input$gene_input!='') {
      matching_ids <- genesToMatchingID[genesToMatchingID$gene_name==input$gene_input,'matching_id']
      selected <- selected & collated$matching_id %in% matching_ids
    }
    
    if (input$variant_input!='') {
      selected <- selected & collated$variant_normalized==input$variant_input
    }
    
    table <- collated[selected,]
    
    if (nrow(table)>0) {
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$data_table <- DT::renderDataTable({
    table <- tableData()
    DT::datatable(table[,c('chemical_normalized','gene_names','variant_id','variant_normalized','chemical_in_pharmgkb','variant_in_pharmgkb','association_in_pharmgkb','in_pharmgkb_paper','citation_count')],
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Chemical','Gene(s)','Variant ID','Variant', 'Chemical in PharmGKB', 'Variant in PharmGKB', 'Association in PharmGKB','In PharmGKB Curated Paper', '# Papers'),
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  dataTableProxy = dataTableProxy('data_table')
  
  output$sentences_table <- DT::renderDataTable({
    if(length(input$data_table_rows_selected)==0) {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    } else {
      table <- tableData()
      row <- table[input$data_table_rows_selected,]
      entries <- sentences[sentences$matching_id==row$matching_id,]
    }
    
    DT::datatable(entries[,c('pubmed_link','pmid_in_pharmgkb', 'journal_short','year', 'subsection' ,'formatted_sentence')],
                  selection = 'none',
                  rownames = FALSE,
                  colnames = c('PMID','PMID in PharmGKB', 'Journal', 'Year', 'Subsection', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  
  chemical_event_val <- reactiveValues(count = 0)
  chemical_piechart_data <- reactive({
    table <- tableData()
    piecounts <- NULL
    if (nrow(table) > 0) {
      piecounts <- aggregate(table$citation_count,by=list(table$chemical_normalized),FUN=sum)
      colnames(piecounts) <- c('label','total_freq')
      piecounts <- piecounts[order(piecounts$total_freq,decreasing=T),]
      
      total <- sum(piecounts$total_freq)
      cutoff <- (1.5/100.0) * total
      piecounts$label <- as.character(piecounts$label)
      piecounts[piecounts$total_freq<cutoff,'label'] <- 'other'
      piecounts$label <- factor(piecounts$label, levels=c('other',as.character(piecounts[piecounts$total_freq>=cutoff,'label'])))
      
      piecounts <- aggregate(piecounts$total_freq,by=list(piecounts$label),FUN=sum)
      colnames(piecounts) <- c('label','total_freq')
    }
    piecounts
  })
  
  output$chemical_piechart <- renderPlotly({
    piecounts <- chemical_piechart_data()
    if (!is.null(piecounts)) {
      sourceName <- paste('chemical_piechart_',chemical_event_val$count,sep='')
      p <- plot_ly(piecounts, labels = ~label, values = ~total_freq, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
        layout(title = paste('Chemicals'),
               showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('chemical_piechart_',chemical_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    piecounts <- chemical_piechart_data()
    if (length(d) > 0 && !is.null(piecounts))
    {
      selected <- piecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "chemical_input",selected=selected)
        chemical_event_val$count <- chemical_event_val$count + 1
      }
    }
  })
  
  
  
  gene_event_val <- reactiveValues(count = 0)
  gene_piechart_data <- reactive({
    table <- tableData()
    table$gene_names <- factor(table$gene_names, c(levels(table$gene_names),'[unknown]'))
    table[table$gene_names=='','gene_names'] <- '[unknown]'
    piecounts <- NULL
    if (nrow(table) > 0) {
      piecounts <- aggregate(table$citation_count,by=list(table$gene_names),FUN=sum)
      if (nrow(piecounts) > 0) {
        colnames(piecounts) <- c('label','total_freq')
        piecounts <- piecounts[order(piecounts$total_freq,decreasing=T),]
        
        total <- sum(piecounts$total_freq)
        cutoff <- (1.5/100.0) * total
        piecounts$label <- as.character(piecounts$label)
        piecounts[piecounts$total_freq<cutoff,'label'] <- 'other'
        piecounts$label <- factor(piecounts$label, levels=c('other',as.character(piecounts[piecounts$total_freq>=cutoff,'label'])))
        
        piecounts <- aggregate(piecounts$total_freq,by=list(piecounts$label),FUN=sum)
        colnames(piecounts) <- c('label','total_freq')
      }
    }
    piecounts
  })
  
  output$gene_piechart <- renderPlotly({
    piecounts <- gene_piechart_data()
    if (!is.null(piecounts)) {
      sourceName <- paste('gene_piechart_',gene_event_val$count,sep='')
      p <- plot_ly(piecounts, labels = ~label, values = ~total_freq, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
        layout(title = paste('Genes'),
               showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('gene_piechart_',gene_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    piecounts <- gene_piechart_data()
    if (length(d) > 0 && !is.null(piecounts))
    {
      selected <- piecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "gene_input", selected=selected)
        gene_event_val$count <- gene_event_val$count + 1
      }
    }
  })
  
  
  variant_event_val <- reactiveValues(count = 0)
  variant_piechart_data <- reactive({
    table <- tableData()
    piecounts <- NULL
    if (nrow(table) > 0) {
      piecounts <- aggregate(table$citation_count,by=list(table$variant_normalized),FUN=sum)
      colnames(piecounts) <- c('label','total_freq')
      piecounts <- piecounts[order(piecounts$total_freq,decreasing=T),]
      
      total <- sum(piecounts$total_freq)
      cutoff <- (1.5/100.0) * total
      piecounts$label <- as.character(piecounts$label)
      piecounts[piecounts$total_freq<cutoff,'label'] <- 'other'
      piecounts$label <- factor(piecounts$label, levels=c('other',as.character(piecounts[piecounts$total_freq>=cutoff,'label'])))
      
      piecounts <- aggregate(piecounts$total_freq,by=list(piecounts$label),FUN=sum)
      colnames(piecounts) <- c('label','total_freq')
    }
    piecounts
  })
  
  output$variant_piechart <- renderPlotly({
    piecounts <- variant_piechart_data()
    if (!is.null(piecounts)) {
      sourceName <- paste('variant_piechart_',variant_event_val$count,sep='')
      p <- plot_ly(piecounts, labels = ~label, values = ~total_freq, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
        layout(title = paste('Variants'),
               showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('variant_piechart_',variant_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    piecounts <- variant_piechart_data()
    if (length(d) > 0 && !is.null(piecounts))
    {
      selected <- piecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "variant_input", selected=selected)
        variant_event_val$count <- variant_event_val$count + 1
      }
    }
  })
  
  
  
  
  
  
  
  output$download_collated_all <- downloadHandler(
    filename = function() {
      return("pgxmine_collated.tsv")
    },
    content = function(file) {
      write.table(collated, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_collated_shown <- downloadHandler(
    filename = function() {
      return("pgxmine_collated_subset.tsv")
    },
    content = function(file) {
      outdata <- tableData()
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_sentences_all <- downloadHandler(
    filename = function() {
      return("pgxmine_sentences.tsv")
    },
    content = function(file) {
      write.table(sentences, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_sentences_shown <- downloadHandler(
    filename = function() {
      return("pgxmine_sentences_selectedbiomarker.tsv")
    },
    content = function(file) {
      if(length(input$data_table_rows_selected)>0) {
        table <- tableData()
        row <- table[input$data_table_rows_selected,]
        entries <- sentences[sentences$matching_id==row$matching_id,]
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$download_sentences_above <- downloadHandler(
    filename = function() {
      return("pgxmine_sentences_multiplebiomarkers.tsv")
    },
    content = function(file) {
      table <- tableData()
      entries <- sentences[sentences$matching_id %in% table$matching_id,]
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  observeEvent(input$chemical_clear, {
    updateSelectizeInput(session, "chemical_input", selected = F)
  })  
  observeEvent(input$gene_clear, {
    updateSelectizeInput(session, "gene_input", selected = F)
  })  
  observeEvent(input$variant_clear, {
    updateSelectizeInput(session, "variant_input", selected = F)
  })  
  
  
}

shinyApp(ui, server)
