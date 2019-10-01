library(bookdown)
library(tinytex)
library(data.table)
library(R.utils)
library(stringr)

projectDir <- "C:/Users/jakel/Documents/R/pgxmine/paper"
setwd(projectDir)
  
mdFile = 'paper.Rmd'
stopifnot(file.exists(mdFile))

if (file.exists('_main.Rmd'))
  file.remove('_main.Rmd')

bookdown::render_book(mdFile, "bookdown::pdf_book", clean = F, encoding = "UTF-8", config_file = "_bookdown.yml")

# Let's check the output is all there
stopifnot(dir.exists('_book'))
stopifnot(file.exists('_book/_main.tex'))
#stopifnot(dir.exists('_bookdown_files/_main_files'))

tempDir <- '_paperover_files'
if (dir.exists(tempDir)) {
  unlink(tempDir,recursive=T)
}
dir.create(tempDir)

#everything <- dir(path=".", pattern="*")
#files <- everything[mapply(everything,FUN=isFile)]
#file.copy(files,tempDir)

images <- dir(path=".", pattern="*.png")
file.copy(images,tempDir)

bibFile <- dir(path=".", pattern="*.bib")
stopifnot(length(bibFile)==1)

bib <- readChar(bibFile, file.info(bibFile)$size)

#bib <- "title={CancerMine: A literature-mined resource for drivers, oncogenes and tumor suppressors in cancer},"
#file.copy(bibFile,paste(tempDir,'/bibliography.bib',sep=''))
bib <- gsub("title=\\{(.*?)\\}","title=\\{\\{\\1\\}\\}",bib)
bib <- gsub("\r","",bib)
bibFileConn <- file(paste(tempDir,'/bibliography.bib',sep=''))
writeLines(bib, bibFileConn)
close(bibFileConn)

templateDir <- 'latex'
#file.copy(templateDir, tempDir, recursive=TRUE)
copyDirectory(templateDir,tempDir)

if (dir.exists('_bookdown_files/_main_files')) {
  copyDirectory('_bookdown_files/_main_files',paste(tempDir,'/_main_files',sep=''))
}

texFile = '_book/_main.tex'
stopifnot(file.exists(texFile))
tex <- readChar(texFile, file.info(texFile)$size)

title <- regmatches(tex, regexpr("\\\\title\\{.*?\\}", tex))
title <- sub("^\\\\title\\{", "", title)
title <- sub("\\}$", "", title)
title <- trim(title)

abstract <- regmatches(tex, regexpr("\\\\section\\{Abstract\\}.*?\\\\section", tex))
abstract <- sub("^\\\\section\\{Abstract\\}", "", abstract)
abstract <- sub("^\\\\label\\{abstract\\}", "", abstract)
abstract <- sub("\\\\section$", "", abstract)
abstract <- trim(abstract)

contents <- regmatches(tex, regexpr("\\\\section\\{Abstract\\}.*?\\\\section\\{References\\}", tex))
if (identical(contents, character(0))) {
  contents <- regmatches(tex, regexpr("\\\\section\\{Abstract\\}.*?\\\\end\\{document\\}", tex))
  contents <- sub("\\\\end\\{document\\}$", "", contents)
} else {
  contents <- sub("\\\\section\\{References\\}$", "", contents)
}
contents <- sub("^\\\\section\\{Abstract\\}.*?\\\\section", "\\\\section", contents)
contents <- trim(contents)



templateFile = paste(tempDir,'/template.tex',sep='')
stopifnot(file.exists(templateFile))
template <- readChar(templateFile, file.info(templateFile)$size)
Encoding(template) <- "UTF-8"

template <- gsub("TITLE",title,template,fixed=T)
template <- gsub("ABSTRACT",abstract,template,fixed=T)
template <- gsub("CONTENTS",contents,template,fixed=T)

template <- gsub("\\citep{","\\cite{",template,fixed=T)
template <- gsub("\\citet{","\\cite{",template,fixed=T)
template <- gsub("\\citep[","\\cite[",template,fixed=T)
template <- gsub("\\citet[","\\cite[",template,fixed=T)


# Clean up awful multi-citation nonsense
# where something like [@percha2018global,@mahmood2016dimex] goes to some weird embedded mix of citep/citet
template <- gsub("\\cite\\[(.*?)\\]\\{(.*?)\\}","\\cite{\\2,\\1}",template)
while (T) {
  oldtemplate <- template
  template <- gsub("\\\\cite{([^}]*)\\\\cite{([^}]*)}","\\\\cite{\\1\\2",oldtemplate,perl=T)
  if (template == oldtemplate) {
    break
  }
}

template <- gsub("\r","",template)

template <- gsub("\\begin{figure}","\\begin{figure}[tb]",template,fixed=T)

template <- gsub("\\\\includegraphics\\{","\\\\includegraphics[width=0.9\\\\textwidth]\\{",template)
template <- gsub("\\\\includegraphics\\[width=.*?\\]\\{","\\\\includegraphics[width=0.95\\\\textwidth]\\{",template)


template <- gsub("\\\\includegraphics\\[(.*?)\\]\\{(.*?)\\}","\\\\begin\\{center\\}\n\\\\includegraphics[\\1]\\{\\2\\}\n\\\\end\\{center\\}",template)


#template <- gsub("α","$ \alpha $",template,fixed=T)
Encoding(template) <- "UTF-8"
template <- gsub("\u03B1","$ \\alpha $",template,fixed=T)
Encoding(template) <- "UTF-8"

outFile <- paste(tempDir,'/main.tex',sep='')
fileConn <- file(outFile)
writeLines(template, fileConn)
close(fileConn)


setwd(tempDir)
#tinytex::pdflatex('main.tex')
pdflatex('main.tex')
#xelatex('main.tex')

setwd(projectDir)


junk <- dir(path=".", pattern="VennDiagram*") 
file.remove(junk)

