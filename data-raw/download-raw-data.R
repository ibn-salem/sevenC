

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

