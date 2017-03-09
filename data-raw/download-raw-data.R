

#-------------------------------------------------------------------------------
# Download an exampel bigWig file from ENCODE
#-------------------------------------------------------------------------------

bwURLstr <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig"
bwTmpFile <- "data-raw/GM12878_Stat1.bigWig"
bwFinalFile <- "inst/extdata/GM12878_Stat1.chr22_1-18000000.bigWig"
dir.create("inst/extdata", recursive = TRUE)

#require(RCurl)
# bwURL <- getURLContent()

download.file(bwURLstr, destfile=bwTmpFile)

bw <- rtracklayer::import(bwTmpFile)
bwSub <- IRanges::subsetByOverlaps(bw, GenomicRanges::GRanges("chr22", IRanges::IRanges(1, 18000000)))

print(object.size(bwSub), unit="Kb")

# save subset of file in inst/extdata as bigWig file
rtracklayer::export(bwSub, con=bwFinalFile, format="bigWig")


#-------------------------------------------------------------------------------
# add CTCF moitif locations in genome
#-------------------------------------------------------------------------------
CTCF_hg19_file <- "data-raw/hg19_MA0139.1.mrkv1.ms"

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

# parse matrix-scan motif:
motifDF <- read.table(CTCF_hg19_file,
                      sep="\t", header=TRUE, comment.char=";")

motif.hg19.CTCF <- GRanges(
  motifDF[,1],
  IRanges(motifDF$start, motifDF$end),
  strand=ifelse(motifDF$strand == "D", "+", "-"),
  motifDF[,c("sequence", "weight", "Pval", "ln_Pval", "sig")],
  seqinfo=hg19seqInfo
)

motif.hg19.CTCF <- sort(motif.hg19.CTCF)

motif.hg19.CTCF <- motif.hg19.CTCF[motif.hg19.CTCF$Pval <= 10^-5]

print(object.size(motif.hg19.CTCF), unit="auto")
print(length(motif.hg19.CTCF))

devtools::use_data(motif.hg19.CTCF, overwrite=TRUE)


#-------------------------------------------------------------------------------
# Download sample loops from Rao et al. 2014 study
#-------------------------------------------------------------------------------
loopURL <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fprimary%2Breplicate%5FHiCCUPS%5Flooplist%5Fwith%5Fmotifs%2Etxt%2Egz"

tmp <- tempfile()

download.file(loopURL, tmp)

df <- read.table(gzfile(tmp), header=TRUE)

# filter for only chromosome 22 in some start region
fltDF <- subset(df, chr1 == 22 &
                  chr2 == 22 &
                  x2 <= 18000000 &
                  y2 <= 18000000)

# save subset of file in inst/extdata as raw .txt file
loopFinalFile <- "inst/extdata/GM12878_HiCCUPS.chr22_1-18000000.loop.txt"

write.table(fltDF, file=loopFinalFile, col.names = TRUE, row.names=FALSE,
            quote=FALSE, sep="\t")

