###################
## Toy example showing R resources to downaload, process and plot genomic features

library(GenomicRanges)
library(FDb.UCSC.snp137common.hg19)

## Example region around Nanog gene
myChr <- "chr12"
myStart <- 7941000
myEnd <- 7949000


## Load genes from UCSC, takes a few mins
txdb <- makeTranscriptDbFromUCSC(genome="hg19", tablename="knownGene")

## Load SNPs from Bioconductor package
snpdb <- FDb.UCSC.snp137common.hg19


## Get SNPs within genes in our region interest
interesting.region <- GRanges(seqnames=myChr, ranges=IRanges(start=myStart, end=myEnd))
interesting.region

tx.ranges <- exons(txdb)
tx.ranges
tx.ranges <- subsetByOverlaps(tx.ranges, interesting.region, type="any", ignore.strand = TRUE)
tx.ranges

snp.ranges <- features(FDb.UCSC.snp137common.hg19) # takes a few mins
snp.ranges
snp.ranges <- subsetByOverlaps(snp.ranges, interesting.region, type="any", ignore.strand = TRUE)
snp.ranges

tx.snp.ranges <- subsetByOverlaps(snp.ranges, tx.ranges, type="any", ignore.strand = TRUE)
tx.snp.ranges

## Make annotation tracks, and plot them
library(Gviz)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = myChr)

txTrack <- GeneRegionTrack(txdb, genome = "hg19", chromosome = myChr,
													 name = "Genes", showId=TRUE, col="black", fill="deepskyblue2", geneSymbol=TRUE)

txSnpTrack <- AnnotationTrack(range=tx.snp.ranges, name="SNPs in transcripts", 
															chromosome=myChr, from=myStart, to=myEnd, genome="hg19", col="red")

snpTrack <- AnnotationTrack(range=snp.ranges, name="all SNPs", col="gray")

plotTracks(list(itrack, gtrack, txTrack, txSnpTrack, snpTrack), from=myStart, to=myEnd)







