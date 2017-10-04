require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(tidyverse)

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
CTCF_hg19_file <- "data-raw/hg19_MA0139.1.mrkv1.ms.gz"

hg19seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

# parse matrix-scan motif:
motifDF <- readr::read_tsv(CTCF_hg19_file, comment = ";")

# filter for p-values <= 10^-6
motifDF <- motifDF[motifDF$Pval <= 10^-6, ]

motif.hg19.CTCF <- GRanges(
  motifDF[[1]],
  IRanges(motifDF$start, motifDF$end),
  strand = ifelse(motifDF$strand == "D", "+", "-"),
  as.data.frame(motifDF[,c("sequence", "weight", "Pval", "ln_Pval", "sig")]),
  seqinfo = hg19seqInfo
)

motif.hg19.CTCF <- sort(motif.hg19.CTCF)


print(object.size(motif.hg19.CTCF), unit = "auto")
print(length(motif.hg19.CTCF))

devtools::use_data(motif.hg19.CTCF, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add subset of CTCF moitif locations on chr22
#-------------------------------------------------------------------------------
motif.hg19.CTCF.chr22 <- motif.hg19.CTCF[seqnames(motif.hg19.CTCF) == "chr22"]
  # & end(motif.hg19.CTCF) <= 18000000]
devtools::use_data(motif.hg19.CTCF.chr22, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add default model parametes trained on best 10 performing TF ChIP-seq data
#-------------------------------------------------------------------------------
model_best_n_file <- "data-raw/v04_screen_TF_lfc.motifSig6_w1000_b1.bestNModelDF.tsv"
modelBest10Avg <- readr::read_tsv(model_best_n_file) %>%
  select(term, estimate_mean) %>%
  rename(estimate = estimate_mean) %>%
  as.data.frame()
devtools::use_data(modelBest10Avg, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add best f1-score cutoff for each TF
#-------------------------------------------------------------------------------
f1ModelDF_file <- "data-raw/v04_screen_TF_lfc.motifSig6_w1000_b1.f1ModelDF.tsv"

cutoffByTF <- readr::read_tsv(f1ModelDF_file) %>%
  select(-max_idx)

devtools::use_data(cutoffByTF, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add default cutoff using optimal f1-score on best 10 performing TF datasets
#-------------------------------------------------------------------------------
topNf1ModelDF_file <- "data-raw/v04_screen_TF_lfc.motifSig6_w1000_b1.topNf1ModelDF.tsv"
cutoffBest10 <- readr::read_tsv(topNf1ModelDF_file)$mean_max_cutoff

devtools::use_data(cutoffBest10, overwrite = TRUE)

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

#-------------------------------------------------------------------------------
# Download sample ChIA-pet data from Tang et al. 2015 Cell
#-------------------------------------------------------------------------------
loopURL <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1872nnn/GSM1872886/suppl/GSM1872886%5FGM12878%5FCTCF%5FPET%5Fclusters.txt.gz"

tmp <- tempfile()
download.file(loopURL, tmp)

df <- read.table(gzfile(tmp), header = FALSE)

# filter for only chromosome 22 in some start region
fltDF <- subset(df, V1 == "chr22" &
                  V4 == "chr22" &
                  V3 <= 18000000 &
                  V6 <= 18000000)

# save subset of file in inst/extdata as raw .txt file
loopFinalFile <- "inst/extdata/ChIA-PET_GM12878_Tang2015.chr22_1-18000000.clusters.txt"

write.table(fltDF, file=loopFinalFile, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")


# #-------------------------------------------------------------------------------
# # Download Capture-Hi-C sample data from Mifsud et al. 2015
# #-------------------------------------------------------------------------------
#
# work_dir <- getwd()
# zipURL <- "http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip"
# zipFile <- "E-MTAB-2323.additional.1.zip"
# inFile <- "TS5_GM12878_promoter-promoter_significant_interactions.txt"
# outFile <- paste0(work_dir, "/inst/extdata/", inFile, ".chr22")
#
# setwd("data-raw")
#
# download.file(zipURL, zipFile)
# system(paste("unzip", zipFile))
# # grep for header and chr22
# system(
#   paste("grep -E -e \"^chr[[:space:]]\" -e \"chr22\" ", inFile, ">", outFile)
#   )
#
# setwd(work_dir)




