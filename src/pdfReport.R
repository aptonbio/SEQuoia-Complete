#' ---
#' title: |
#'   | \vspace{2cm} \LARGE{SEQuoia Analysis Report}
#' header-includes:
#' - \usepackage{titling}
#' - \usepackage{fontspec}
#' - \setmainfont{FreeSans}
#' - \usepackage{booktabs}
#' - \usepackage[table]{xcolor}
#' - \pretitle{\vspace{5cm}\begin{center}\LARGE\includegraphics[width=12cm]{/opt/biorad/src/vendor-logo.png}\\[\bigskipamount]}
#' - \posttitle{\end{center}\newpage}
#' output: 
#'   pdf_document:
#'     latex_engine: xelatex
#'     toc: true
#' always_allow_html: yes
#' ---


#' \centering
#' \raggedright
#' \newpage

#+ setup, include=FALSE

#basics
library(knitr)
library(kableExtra)
library(dplyr) 
library(data.table)
library(ggplot2)
library(tibble)
library(plotly)
library(fastqcr)
library(rlist)

#muting warnings
options(warn=-1)

#setwd("/")
#dictate kable styling options
kableStyle <- c("striped", "condensed", "hover", "responsive")

#fastqc
fastqcDir <- paste(base_dir, "fastqc", sep="/")
fastqcDirExists <- dir.exists(fastqcDir)
write(paste("fastqcDirExists: ", fastqcDirExists), stderr())

#debarcoding
debarcodeDir <- paste(base_dir, "debarcode", sep="/")
debarcodeDirExists <- dir.exists(debarcodeDir)
write(paste("debarcodeDirExists: ", debarcodeDirExists), stderr())

#trimming
trimDir <- paste(base_dir, "cutadapt", sep="/")
trimDirExists <- dir.exists(trimDir)
write(paste("trimDirExists: ", trimDirExists), stderr())

#alignments
alignmentDir <- paste(base_dir, "alignment", sep="/")
alignmentDirExists <- dir.exists(alignmentDir)
write(paste("alignmentDirExists: ", alignmentDirExists), stderr())

#deduplication
dedupDir <- paste(base_dir, "dedup", sep="/")
dedupDirExists <- dir.exists(dedupDir)
write(paste("dedupDirExists: ", dedupDirExists), stderr())

#counts
countsDir <- paste(base_dir, "counts", sep="/")
countsDirExists <- dir.exists(countsDir)
write(paste("countsDirExists: ", countsDirExists), stderr())

#' `r if(fastqcDirExists) { "# Read QC" }`
#+ eval=fastqcDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F

parseQcFile <- function(fileName)
{
  qc <- qc_read(fileName, modules="Per base sequence quality")
  return(qc$per_base_sequence_quality)
}

qcFiles <- list.files(fastqcDir, full.names=TRUE)[grepl("zip",list.files(fastqcDir))]

r1 <- qcFiles[grepl("_R1_",qcFiles)]
r2 <- qcFiles[grepl("_R2_",qcFiles)]

if (length(r1) == 0 || length(r2) == 0) {
  stop("cannot find QC files containing R1/R2 in file name")
}

r1Qc <- sapply(r1, parseQcFile, simplify=F)

r1QcDf <- list.cbind(r1Qc)

p1 <- plot_ly(width = 700) %>% 
  layout(
    yaxis = list(title = "Q Score", range = c(0, 40)),
    xaxis = list(title = "Position in read (bp)",
                 tickvals = seq(1, dim(r1QcDf)[1], by=2),
                 tickmode = "array",
                 ticktext = r1Qc[[1]]$Base[seq(1, dim(r1QcDf)[1], by=2)],
                 tickangle = 90,
                 ticks = "outside"))

for(i in grep("Mean",names(r1QcDf)))
{
  p1 <- add_trace(p1, x = ~c(1:(length(r1QcDf[,1]))), y = r1QcDf[,i], type = 'scatter', mode ='lines', name = gsub("_fastqc.zip.Mean", " ", basename(names(r1QcDf[i]))))
}

p1 <- p1 %>% layout(
  shapes = list(
    list(type = "rect",
         fillcolor = "red",
         line = list(color="red"),
         opacity = 0.3,
         x0=0, x1=length(r1QcDf[,1]), xref="x",
         y0=0, y1=20, yref="y",
         layer = "below"),
    list(type = "rect",
         fillcolor = "yellow",
         line = list(color="yellow"),
         opacity = 0.3,
         x0=0, x1=length(r1QcDf[,1]), xref="x",
         y0=20, y1=28, yref="y",
         layer = "below"),
    list(type = "rect",
         fillcolor = "green",
         line = list(color="green"),
         opacity = 0.3,
         x0=0, x1=length(r1QcDf[,1]), xref="x",
         y0=28, y1=40, yref="y",
         layer = "below")
    ),
    legend = list(
      orientation = 'h',
      x = 10,
      y = -0.15)
)

if (length(r2) > 0) {

  r2Qc <- sapply(r2, parseQcFile, simplify=F)
  r2QcDf <- list.cbind(r2Qc)
  p2 <- plot_ly(width = 700) %>% 
    layout(
      yaxis = list(title = "Q Score", range = c(0, 40)),
      xaxis = list(title = "Position in read (bp)",
       tickvals = seq(1, dim(r2QcDf)[1], by=2),
       tickmode = "array",
       ticktext = r1Qc[[1]]$Base[seq(1, dim(r2QcDf)[1], by=2)],
       tickangle = 90,
       ticks = "outside"))

  for(i in grep("Mean", names(r2QcDf)))
  {
    p2 <- add_trace(p2, x = ~c(1:(length(r2QcDf[,1]))), y = r2QcDf[,i], type = 'scatter', mode ='lines', name = gsub("_fastqc.zip.Mean", " ", basename(names(r2QcDf[i]))))
  }

  p2 <- p2 %>% layout(
    shapes = list(
      list(type = "rect",
     fillcolor = "red",
     line = list(color="red"),
     opacity = 0.3,
     x0=0, x1=length(r2QcDf[,1]), xref="x",
     y0=0, y1=20, yref="y",
     layer = "below"),
      list(type = "rect",
     fillcolor = "yellow",
     line = list(color="yellow"),
     opacity = 0.3,
     x0=0, x1=length(r2QcDf[,1]), xref="x",
     y0=20, y1=28, yref="y",
     layer = "below"),
      list(type = "rect",
     fillcolor = "green",
     line = list(color="green"),
     opacity = 0.3,
     x0=0, x1=length(r2QcDf[,1]), xref="x",
     y0=28, y1=40, yref="y",
     layer = "below")
    ),
    legend = list(
      orientation = 'h',
      x = 10,
      y = -0.15)
  )
}

#' `r if(exists("p1")) { "## Read 1" }`
#+ eval=exists("p1"), echo=FALSE, fig.asp=1, fig.align="center", message=F
p1

#' `r if(exists("p2")) { "## Read 2" }`
#+ eval=exists("p2"), echo=FALSE, fig.asp=1, fig.align="center", message=F
p2

#' \newpage

#' `r if(debarcodeDirExists) { "# UMI Parsing" }`
#+ eval=debarcodeDirExists, echo=FALSE, fig.width=8, fig.height=5, fig.align="left"

deb <- read.table(list.files(debarcodeDir, full.names=TRUE)[0], fill=T)
inputReads <- as.numeric(as.character(deb$V3[1]))
validBcReads <- as.numeric(as.character(deb$V3[2]))
invalidBcReads <- inputReads-validBcReads

#create data frame
df <- data.frame(
  Metric = c("Input Reads", "Reads with Valid UMI", "% Reads with Valid UMI"),
  Value = prettyNum(c(inputReads, validBcReads, signif(validBcReads/inputReads, 3) * 100), big.mark = ",", scientific = F),
  stringsAsFactors = FALSE
)

#create kable output
kable(df, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#plotly
dfx <- data.frame(" "="barcode", "valid"=validBcReads, "invalid"=inputReads-validBcReads)
pl <- plot_ly(dfx, x = ~valid, y = ~" ", type = "bar", name="Valid", orientation = "h", hoverinfo = 'text', text = ~paste("Valid Reads: ", valid)) %>% 
  add_trace(x = ~invalid, name = "Invalid", hoverinfo = 'text', text = ~paste("Invalid Reads: ", invalid)) %>% 
  layout(barmode = 'stack', xaxis = list(title = "Reads"), yaxis = list(title = ""))

#create plot
pl

#' \newpage

#' `r if(trimDirExists) { "# Read Trimming" }`
#+ eval=trimDirExists, echo=FALSE, fig.asp=0.75, fig.align="center"

#Import metadata from cutadapt and format it
rt <- read.table(list.files(trimDir, full.names=TRUE)[0], skip=7, nrows=7, fill=T, sep=":") 
names(rt) <- c("Metric","Value")
rt$Value <- as.numeric(gsub(",","",unlist(lapply(strsplit(as.character(rt$Value), split="\\s+"), `[[`, 2)))) #this is gross, i'm sorry for nesting 6 functions

#generate table for plotting
dfx <- data.frame(" "="Reads",
                  "Reads Input" = rt$Value[rt$Metric=="Total reads processed"],
                  "Reads Too Short" = rt$Value[rt$Metric=="Reads that were too short"],
                  "Reads Written" = rt$Value[rt$Metric=="Reads written (passing filters)"], check.names=FALSE)

pl <- plot_ly(dfx, x = ~`Reads Too Short`, y = ~" ", type = "bar", name = "Reads Too Short", orientation = "h", hoverinfo = "text", text = ~paste("Reads Too Short: ", dfx$`Reads Too Short`)) %>%
  add_trace(x = ~`Reads Written`, name = "Reads Written", hoverinfo = "text", text = ~paste("Reads Written: ", dfx$`Reads Written`)) %>%
  layout(barmode = "stack", xaxis = list(title = "Reads"), yaxis = list(title = ""), legend = list(orientation = 'h', x=0.3, y=-0.4))

#write table
rt$Value <- prettyNum(rt$Value, big.mark = ",", scientific = F) #this is gross, i'm sorry for nesting 6 functions
kable(rt, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


#' \newpage

#' `r if(alignmentDirExists) { "# Alignment" }`
#+ eval=alignmentDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F

#Using a subset of PICARD outputs
alignmentFiles <- list.files(alignmentDir, full.names=TRUE)
rt <- read.table(alignmentFiles[grepl("rna_metrics.txt", alignmentFiles)], nrows=1, header=T, fill=T)
starReport <- read.table(alignmentFiles[grepl("Log.final.out", alignmentFiles)], sep="|", fill=T)
starReport$V2 <- gsub("\t","",starReport$V2)

#parse out relevant stuff from STAR report
inputReads <- as.numeric(starReport$V2[grep("Number of input reads",starReport$V1)])
uniqueMapped <- as.numeric(starReport$V2[grep("Uniquely mapped reads number",starReport$V1)])
multiMapped <- as.numeric(starReport$V2[grep("Number of reads mapped to multiple loci",starReport$V1)])
tooManyMapped <- as.numeric(starReport$V2[grep("Number of reads mapped to too many loci",starReport$V1)])
unmapped <- inputReads-uniqueMapped-multiMapped

df <- data.frame("Reads Input" = inputReads,
                 "Uniquely Mapped Reads" = uniqueMapped,
                 "Multi-mapped Reads" = multiMapped,
                 "Reads mapped to too many loci" = tooManyMapped,
                 "Unmapped Reads" = unmapped,
                 "PF Bases" = rt$PF_BASES,
                 "PF Aligned Bases" = rt$PF_ALIGNED_BASES,
                 "Coding Bases" = rt$CODING_BASES,
                 "UTR Bases" = rt$UTR_BASES,
                 "Intronic Bases" = rt$INTRONIC_BASES,
                 "Intergenic Bases" = rt$INTERGENIC_BASES,
                 "Ribosomal Bases" = rt$RIBOSOMAL_BASES,
                 "Median CV Coverage" = rt$MEDIAN_CV_COVERAGE,
                 "Median 5' Bias" = rt$MEDIAN_5PRIME_BIAS,
                 "Median 3' Bias" = rt$MEDIAN_3PRIME_BIAS,
                 "Median 5' to 3' Bias" = rt$MEDIAN_5PRIME_TO_3PRIME_BIAS,
                 "% Stranded" = rt$PCT_CORRECT_STRAND_READS * 100,
		 "% rRNA bases" = (rt$RIBOSOMAL_BASES/rt$PF_BASES)*100,
                 check.names= F)
df <- as.data.frame(t(df)) %>% rownames_to_column()
names(df) <- c("Metric", "Value")

dfx <- data.frame(" "="Bases",
                  "Coding" = df$Value[which(df$Metric=="Coding Bases")],
                  "UTR" = df$Value[which(df$Metric=="UTR Bases")],
                  "Intronic" = df$Value[which(df$Metric=="Intronic Bases")],
                  "Intergenic" = df$Value[which(df$Metric=="Intergenic Bases")],
                  "Ribosomal" = df$Value[which(df$Metric=="Ribosomal Bases")], check.names=F)

pl <- plot_ly(dfx, x = ~`Coding`, y= ~" ", type = "bar", name="Coding", orientation = "h", hoverinfo = 'text', text = ~paste("Coding: ", dfx$`Coding`)) %>% 
  add_trace(x = ~`UTR`, name = "UTR", hoverinfo = 'text', text = ~paste("UTR: ", round(dfx$`UTR`))) %>%
  add_trace(x = ~`Intronic`, name = "Intronic", hoverinfo = 'text', text = ~paste("Intronic: ", round(dfx$`Intronic`))) %>%
  add_trace(x = ~`Intergenic`, name = "Intergenic", hoverinfo = 'text', text = ~paste("Intergenic: ", round(dfx$`Intergenic`))) %>%
  add_trace(x = ~`Ribosomal`, name = "Ribosomal", hoverinfo = 'text', text = ~paste("Ribosomal: ", round(dfx$`Ribosomal`))) %>%
  layout(barmode = 'stack', xaxis = list(title = "Aligned Bases"), yaxis = list(title = ""), legend = list(orientation = 'h', x=0.3, y=-0.4))

#create table
df$Value <- prettyNum(df$Value, big.mark = ",", scientific = F)
kable(df, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

rt_cov <- read.table(alignmentFiles[grepl("rna_metrics.txt", alignmentFiles)], skip = 10, header=T, fill=T)

cov <- plot_ly(width = 700) %>% 
  layout(
    yaxis = list(title = "Normalized Coverage", range = c(0, max(rt_cov$All_Reads.normalized_coverage+0.1*rt_cov$All_Reads.normalized_coverage))),
    xaxis = list(title = "Normalized Position",
                 tickvals = seq(0, 100, by=2),
                 tickmode = "array",
                 ticktext = seq(0, 100, by=2),
                 tickangle = 90,
                 ticks = "outside"))

cov <- add_trace(cov, x = ~c(1:101), y = rt_cov$All_Reads.normalized_coverage, type = 'scatter', mode ='lines')

#' \newpage

#' `r if(exists("cov")) { "## Transcript Coverage" }`
#+ eval=exists("cov"), echo=FALSE, fig.asp=1, fig.align="center", message=F
cov

#' \newpage

#' `r if(exists("pl")) { "## Base Distribution" }`
#+ eval=exists("pl"), echo=FALSE, fig.asp=1, fig.align="center", message=F
pl

#' \newpage

#' `r if(dedupDirExists) { "# Deduplication" }`
#+ eval=dedupDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
dedupLog = list.files(dedupDir, full.names=TRUE)[0]
umisObserved <- as.numeric(system(paste('grep -F "#umis"', dedupLog, "| cut -d' ' -f5"), intern=T))
inputAlignments <- as.numeric(system(paste('grep "Input Reads"', dedupLog, "| cut -d' ' -f7"), intern=T))
outputAlignments <- as.numeric(system(paste('grep "reads out"', dedupLog, "| cut -d: -f4"), intern=T))
meanUmiPerPos <- as.numeric(system(paste('grep "Mean number of unique UMIs per position"', dedupLog, "| cut -d: -f4"), intern=T))
maxUmiPerPos <- as.numeric(system(paste('grep "Max. number of unique UMIs per position"', dedupLog, "| cut -d: -f4"), intern=T))
uniqInputReads <- as.numeric(system(paste('grep "unique_input_reads"', dedupLog, "| cut -d ' ' -f2"), intern=T))
uniqOutputReads <- as.numeric(system(paste('grep "unique_output_reads"', dedupLog, "| cut -d ' ' -f2"), intern=T))

df <- data.frame("Total input alignments" = inputAlignments,
                 "Total output alignments" = outputAlignments,
                 "Unique UMIs observed" = umisObserved,
                 "Average UMIs per position" = meanUmiPerPos,
                 "Maximum UMIs per position" = maxUmiPerPos,
                 "Unique Input Reads" = uniqInputReads,
                 "Unique Output Reads" = uniqOutputReads,
                 "% PCR Duplicates" = (1 - (uniqOutputReads / uniqInputReads)) * 100,
                 check.names= F)

df <- as.data.frame(t(df)) %>% rownames_to_column()
names(df) <- c("Metric", "Value")

df$Value <- prettyNum(df$Value, big.mark = ",", scientific=FALSE)
kable(df, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage

#' `r if(countsDirExists) { "# Transcriptome" }`
#+ eval=countsDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
#parse genecounts summary for longRNA

countsFiles <- list.files(countsDir, full.names=TRUE)

longRNAcounts <- read.table(countsFiles[grepl("longRNAs.*[:.:]summary", countsFiles)[0]], skip=1)
colnames(longRNAcounts) <- c("Result", "Count")
longRNAcounts$Result <- gsub("_", " ", longRNAcounts$Result)
longRNAcounts <- rbind(data.frame(Result="Total Alignments", Count=sum(longRNAcounts$Count)), longRNAcounts)

#parse genecounts summary for miRNA
miRNAcounts <- read.table(countsFiles[grepl("miRNAs.*[:.:]summary", countsFiles)[0]], skip=1)
colnames(miRNAcounts) <- c("Result", "Count")
miRNAcounts$Result <- gsub("_", " ", miRNAcounts$Result)
miRNAcounts <- rbind(data.frame(Result="Total Alignments", Count=sum(miRNAcounts$Count)), miRNAcounts)

#handle biotypes
countLong <- read.table(countsFiles[grepl("longRNAs.*[:.:]gene_counts$", countsFiles)[0]], sep="\t", header=T, col.names=c("Gene", "Chr", "Start", "End", "Strand", "Length", "Count"))
countMicro <- read.table(countsFiles[grepl("miRNAs.*[:.:]gene_counts$", countsFiles)[0]], sep="\t", header=T, col.names=c("Gene", "Chr", "Start", "End", "Strand", "Length", "Count"))
countMicro$gene_biotype <- "miRNA"
biotypes <- read.table(biotypes_file, sep="\t", header=T)

countLong <- left_join(countLong, biotypes, by = c("Gene" = "gene_id"))
countAll <- rbind(countLong, countMicro)
countByBiotype <- countAll %>% filter(!is.na(gene_biotype)) %>% group_by(gene_biotype) %>% summarise(count = sum(Count)) %>% arrange(-count)
countByBiotype$gene_biotype <- factor(countByBiotype$gene_biotype, levels = unique(countByBiotype$gene_biotype)[order(countByBiotype$count, decreasing = TRUE)])
countByBiotype <- countByBiotype %>% filter(count > 0) #filter biotypes with no counts

#create plot with labels above bars; plotly handles autoscaling
pl <- plot_ly(countByBiotype,
              x=~gene_biotype,
              y=~count,
              text=~prettyNum(count, big.mark = ",", scientific=FALSE),
              textposition='outside',
              type='bar')

pl <- plot_ly(countByBiotype, x=~gene_biotype, y=~count, type='bar', width=700)

countByBiotype$count <- prettyNum(countByBiotype$count, big.mark = ",", scientific=FALSE)

#render tables and plots on separate pages
#' `r if(exists("longRNAcounts")) { "## longRNA Counts" }`
#+ eval=exists("longRNAcounts"), echo=FALSE, fig.asp=1, fig.align="center", message=F
kable(longRNAcounts, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage

#' `r if(exists("miRNAcounts")) { "## miRNA Counts" }`
#+ eval=exists("miRNAcounts"), echo=FALSE, fig.asp=1, fig.align="center", message=F
kable(miRNAcounts, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage

#' `r if(exists("countByBiotype")) { "## Gene Biotypes" }`
#+ eval=exists("countByBiotype"), echo=FALSE, fig.asp=1, fig.align="center", message=F
kable(countByBiotype, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#render plot
pl

#' \newpage

#' `r if(TRUE) { "# Pipeline Metadata" }`
#+ eval=TRUE, echo=FALSE, fig.asp=1, fig.align="center", message=F, results="asis", warn=F
env <- Sys.getenv(c("FASTQC_VERSION","STAR_VERSION","BEDTOOLS_VERSION","PICARD_VERSION","UMI_TOOLS_VERSION","SUBREAD_VERSION","SAMBAMBA_VERSION"))
env <- as.data.frame(env, stringsAsFactors=FALSE) %>% tibble::rownames_to_column()
umi_tools_version <- system("umi_tools --version", intern=T)
umi_tools_version <- strsplit(umi_tools_version, ":")[[1]][2]
env[which(env$rowname=="UMI_TOOLS_VERSION"), 2] = gsub(" ", "", umi_tools_version)
anno_version <- read.table(anno_file, comment.char="", fill=T, sep=",")
anno_source <- gsub("#!annotation-source ", "", anno_version$V1[grep("annotation-source", anno_version$V1)])
localVars <- data.frame(rowname = c("Reference Genome", "Annotation Source", "UMI Aware", "ERCC"), env = c(referenceGenome, anno_source, dedupDirExists, isErcc), stringsAsFactors=FALSE)
env <- rbind(localVars, env)
env[nrow(env) + 1,] = list("Report Generated", paste(as.character(Sys.time()), "UTC"))
colnames(env) <- NULL
kable(env, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))
