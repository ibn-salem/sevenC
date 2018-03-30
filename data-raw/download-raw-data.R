require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(tidyverse)

#-------------------------------------------------------------------------------
# Download an exampel bigWig file from ENCODE
#-------------------------------------------------------------------------------

bwURLstr <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig"
bwTmpFile <- "data-raw/GM12878_Stat1.bigWig"
bwFinalFile <- "inst/extdata/GM12878_Stat1.chr22_1-30000000.bigWig"
dir.create("inst/extdata", recursive = TRUE)

#require(RCurl)
# bwURL <- getURLContent()

download.file(bwURLstr, destfile = bwTmpFile)

# read bigWig file as GRanges
bw <- rtracklayer::import(bwTmpFile)

# get subset on chr22 1:30000000
bwSub <- IRanges::subsetByOverlaps(
  bw, GenomicRanges::GRanges("chr22", IRanges::IRanges(1, 30000000))
)

print(object.size(bwSub), unit = "Mb")

# save subset of file in inst/extdata as bigWig file
rtracklayer::export(bwSub, con = bwFinalFile, format = "bigWig")
# rtracklayer::export(bwChr22, con = bwFinalFile, format = "bigWig")

#*******************************************************************************
# add CTCF moitif locations in human genome from JASPAR motif tracks      ------
#*******************************************************************************

# define link to TSV file from JASPAR tracks
JASPAR_HG19_CTCF_URL <- "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg19/tsv/MA0139.1.tsv.gz"

# define minimal significance threshold as p-value
MOTIF_PVAL <- 2.5 * 1e-06

# Parse CTCF motif sites from JASPAR track

# define colnames and parse TSV file from JASPAR
# header: chr `start (1-based)`   end `rel_score * 1000` `-1 * log10(p_value) * 100` strand
col_names = c("chr", "start", "end", "name", "score", "log10_pval_times_100", "strand")
allMotifDF <- read_tsv(JASPAR_HG19_CTCF_URL, col_names = col_names, skip = 1,
                       col_type = cols(
                         chr = col_character(),
                         start = col_integer(),
                         end = col_integer(),
                         name = col_character(),
                         score = col_integer(),
                         log10_pval_times_100 = col_integer(),
                         strand = col_character()
                       )) %>%
  mutate(log10_pval = log10_pval_times_100 / 100)

motifDF <- allMotifDF %>%
  filter(log10_pval >= -log10(MOTIF_PVAL))

# take seqinfo object from hg19 TxDb package
seqInfoHg19 <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)


# build GRanges object
motif.hg19.CTCF <- GRanges(motifDF$chr, IRanges(motifDF$start, motifDF$end),
                   strand = motifDF$strand,
                   score = motifDF$log10_pval,
                   seqinfo = seqInfoHg19)
# sort
motif.hg19.CTCF <- sort(motif.hg19.CTCF)

# print size and number
print(object.size(motif.hg19.CTCF), unit = "auto")
print(length(motif.hg19.CTCF))

devtools::use_data(motif.hg19.CTCF, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add subset of CTCF moitif locations on chr22
#-------------------------------------------------------------------------------
motif.hg19.CTCF.chr22 <- motif.hg19.CTCF[seqnames(motif.hg19.CTCF) == "chr22"]

devtools::use_data(motif.hg19.CTCF.chr22, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add example ChIP-seq coverage
#-------------------------------------------------------------------------------
motif.hg19.CTCF.chr22.cov <- sevenC::addCovToGR(
  motif.hg19.CTCF.chr22,
  bwFinalFile)
devtools::use_data(motif.hg19.CTCF.chr22.cov, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add default model parametes trained on best 10 performing TF ChIP-seq data
#-------------------------------------------------------------------------------
model_best_n_file <- "data-raw/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1.bestNModelDF.tsv"

modelBest10Avg <- readr::read_tsv(model_best_n_file) %>%
  select(term, estimate_mean) %>%
  rename(estimate = estimate_mean) %>%
  as.data.frame()
devtools::use_data(modelBest10Avg, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add TF specific models as data frame
#-------------------------------------------------------------------------------
TFspecific_ModelDF_meta_file = "data-raw/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1.TFspecific_ModelDF_meta.tsv"
TFspecificModels <- read_tsv(TFspecific_ModelDF_meta_file)
devtools::use_data(TFspecificModels, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add best f1-score cutoff for each TF
#-------------------------------------------------------------------------------
f1ModelDF_file <- "data-raw/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1.f1ModelDF.tsv"

cutoffByTF <- readr::read_tsv(f1ModelDF_file) %>%
  select(-max_idx)

devtools::use_data(cutoffByTF, overwrite = TRUE)

#-------------------------------------------------------------------------------
# add default cutoff using optimal f1-score on best 10 performing TF datasets
#-------------------------------------------------------------------------------
topNf1ModelDF_file <- "data-raw/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1.topNf1ModelDF.tsv"
cutoffBest10 <- readr::read_tsv(topNf1ModelDF_file)$mean_max_cutoff

devtools::use_data(cutoffBest10, overwrite = TRUE)

#-------------------------------------------------------------------------------
# Download sample loops from Rao et al. 2014 study
#-------------------------------------------------------------------------------
loopURL <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fprimary%2Breplicate%5FHiCCUPS%5Flooplist%5Fwith%5Fmotifs%2Etxt%2Egz"

tmp <- tempfile()

download.file(loopURL, tmp)

df <- read.table(gzfile(tmp), header = TRUE)

# filter for only chromosome 22 in some start region
fltDF <- subset(df, chr1 == 22 &
                  chr2 == 22 &
                  x2 <= 30000000 &
                  y2 <= 30000000)

# save subset of file in inst/extdata as raw .txt file
loopFinalFile <- "inst/extdata/GM12878_HiCCUPS.chr22_1-30000000.loop.txt"

write.table(fltDF, file = loopFinalFile, col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")

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
                  V3 <= 30000000 &
                  V6 <= 30000000)

# save subset of file in inst/extdata as raw .txt file
loopFinalFile <- "inst/extdata/ChIA-PET_GM12878_Tang2015.chr22_1-30000000.clusters.txt"

write.table(fltDF, file = loopFinalFile, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")


