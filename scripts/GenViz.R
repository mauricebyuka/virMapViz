#! /usr/bin/Rscript

### Load required libraries
suppressMessages(library(data.table))
suppressMessages(library(Gviz))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(Rsamtools))
suppressMessages(library(tools))

### Set variables
args = commandArgs(trailingOnly=TRUE)
dataValues = args[1]

print('reading data')
data <- read.csv(dataValues, header = TRUE, sep = '\t')

sampleName <- data$sample
out <- data$out_dir
virus <- data$virus
refFile <- data$re_fasta
bamFile <- data$bam
bedFile <- data$bedfile
gffFile <- data$gff # must be GFF3 file from NCBI GenBank/RefSeq (or equivalent)
myChr <- data$chromosome # must match column 1 of gffFile and bedFile
seqEnd <- data$seqEnd

# Coordinates of areas to be plotted
full_Start = 1
full_End = seqEnd
zoomStart <- data$zoomStart
zoomEnd <- data$zoomEnd

graph_dir <- file.path(out, 'graphs')
if (!file.exists(graph_dir)) {
  dir.create(graph_dir)
}

graph1 = paste0(sampleName, '_full_coverage.png')
# graph2 = paste0(sampleName, '_zoom_', zoomStart, '-',zoomEnd, '_coverage.png')
# graph3 = paste0(sampleName, '_full_alignment.png')
# graph4 = paste0(sampleName, '_zoom_', zoomStart, '-',zoomEnd, '_alignment.png')

# PLOTTING

##Add genome axis
genomeTrack <- GenomeAxisTrack(col="black") 

### read in Feature annotation file and set up GeneRegionTrack

fileExt = file_ext(gffFile)
gffExts <- c('gff', 'gff1', 'gff2', 'gff3')

if (fileExt %in% gffExts) {
  virusRange <- makeTxDbFromGFF(gffFile, format  = fileExt)  
} else {
  virusRange <- fread(gffFile, col.names = c('chromosome', 'start', 'end', 'symbol'))
}


print('Setting gene track')


genesTrack <- GeneRegionTrack(virusRange, genome = virus, chromosome = myChr, name = "Genes", col="black", fill="green", stacking="pack", shape="smallArrow", background.title = "darkblue", background.panel = "#FFFEDB", options(ucscChromosomeNames=FALSE), transcriptAnnotation = "symbol", just.group="above", fontcolor.group="black", fontsize.group=15, shape="arrow") #squish #dense

### read in bedgraph dataset and set up DataTrack

print('Reading bed')

file1 <- fread(bedFile, col.names = c('chromosome', 'start', 'end', 'value'))
coverageTrack <- DataTrack(range = file1, chromosome=myChr, genome = virus, fill = "#46a2da", col = "blue", options(ucscChromosomeNames=FALSE), col.axis="black", background.title = "transparent")#, ylim=c(0,max1))

sTrack <- SequenceTrack(refFile)

###Full-size coverage plot ###
fullCovPng <-  file.path(out, 'graphs', graph1)
if (!file.exists(fullCovPng)) {
  png(fullCovPng , width = 1200)
  plotTracks(list(genomeTrack, genesTrack, coverageTrack, sTrack), from = full_Start , to = full_End, sizes=c(0.08, 0.12, 0.20, 0.50), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=0.9, collapse=FALSE)
  dev.off()
}



### Coverage plot for zoomed area
# zoomCovPng <-  file.path(out, 'graphs', graph2)
# if (!file.exists(zoomCovPng)) {
#   png(zoomCovPng, width = 1200)
#   plotTracks(list(genomeTrack, genesTrack, coverageTrack), from = zoomStart, to = zoomEnd, sizes=c(0.06, 0.06, 0.30), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=0.9, collapse=FALSE)
#   dev.off()
# }


### GENERATE ALIGNMENT PLOTS ###

### Load the sequence track & index bam file & Load the alignment track
# alignName = paste0('Alignment:', myChr, '-', sampleName)
# seqTrack <- SequenceTrack(refFile)
# indexBam(bamFile)
# alignTrack <- AlignmentsTrack(bamFile, isPaired = TRUE, name = alignName)
# 

### Full-size alignment plot
# fullAlignPng <-  file.path(out, 'graphs', graph3)
# if (!file.exists(fullAlignPng)){
#   png(fullAlignPng, width = 1200)
#   plotTracks(c(genomeTrack, genesTrack, seqTrack, alignTrack), chromosome = myChr, from = full_Start , to = full_End, options(ucscChromosomeNames=FALSE))
#   dev.off()
# }

### Zoom alignment plot
# zoomAlignPng <-  file.path(out, 'graphs', graph4)
# if (!file.exists(zoomAlignPng)) {
#   png(zoomAlignPng, width = 1200)
#   plotTracks(c(genomeTrack, genesTrack, seqTrack, alignTrack), chromosome = myChr, from = zoomStart, to = zoomEnd, options(ucscChromosomeNames=FALSE), fill.coverage="#46a2da", col.deletion="#FF0000", col.gap="#FFA500", background.title = "darkblue", sizes=c(0.06, 0.08, 0.05, 0.30))
#   dev.off()
# }



