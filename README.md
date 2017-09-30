# chromloop
R package to predict chromatin looping interactions from ChIP-seq data and 
sequnece motifs.

Folding of genomic DNA is an importent feature of
  eukaryotic genomes and can bring regulatory sequences, such as enhancers or 
  transcription factor binding sites, in close physical proximity of regulated 
  target genes.
  Here, we provide a tool that uses protein binding signals from ChIP-seq and
  sequence motif information to predict chromatin looping events. Cross-linking 
  of proteins binding close to loop anchors result in ChIP-seq signals at both 
  anchor loci. These signals are used at CTCF motif pairs together with their 
  distance and orientation to each other to predict whether they interact or not.
  The resulting chromatin loops can be used to associate enhancers or 
  transcription factor binding sites in form of ChIP-seq peaks to regulated 
  target genes.
  
## Intallation

The chromloop package depends on the follwing R packages from Bioconductor:

- rtracklayer (>= 1.34.1),
- InteractionSet (>= 1.2.0),

```R
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer", "InteractionSet")
```

Install chromloop from github using devtools:

```R
#install.packages("devtools")
devtools::install_github("ibn-salem/chromloop")
```

## Basic usage example
Here we show how to use the package to predict chromatin looping interactions 
among CTCF moif locations on chromosome 22. 
As input only a single bigWig file is used from a STAT1 ChIP-seq experiment 
in human GM12878 cells. 

### Get motif pairs
```R
library(InteractionSet)
library(chromloop)

# load provided CTCF motifs
motifs <- motif.hg19.CTCF.chr22

# get motifs pairs
gi <- prepareCandidates(motifs, maxDist = 10^6, scoreColname = "sig")
```

### Add ChIP-seq data and compute correaltion
```R

# use example ChIP-seq bigWig file
bigWigFile <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig", 
  package = "chromloop")

# add ChIP-seq coverage and compute correaltion at motif pairs
gi <- addCor(gi, bigWigFile)
```

###  Predict loops

```R
loops <- predLoops(gi)
```

For more detailed usage instructions, see the package 
[vignette](vignettes/chromloop.Rmd).

