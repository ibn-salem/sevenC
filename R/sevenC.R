#' Computational chromosome conformation capture by correlation of ChIP-seq at
#' CTCF motifs (7C)
#'
#' Chromatin looping is an essential feature of eukaryotic genomes and
#' can bring regulatory sequences, such as enhancers or transcription factor
#' binding sites, in the close physical proximity of regulated target genes.
#' Here, we provide sevenC, an R package that uses protein binding signals from
#' ChIP-seq and sequence motif information to predict chromatin looping events.
#' Cross-linking of proteins that bind close to loop anchors result in ChIP-seq
#' signals at both anchor loci. These signals are used at CTCF  motif pairs
#' together with their distance and orientation to each other to predict
#' whether they interact or not.
#' The resulting chromatin loops might be used to associate enhancers or
#' transcription factor binding sites (e.g., ChIP-seq peaks) to regulated
#' target genes.
#'
#' @name sevenC
"_PACKAGE"
