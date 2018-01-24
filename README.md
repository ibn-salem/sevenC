# sevenC
R package to predict chromatin looping interactions from ChIP-seq by computational 
chromosome conformation capture by correlation of ChIP-seq at CTCF motifs (7C).

Chromatin looping is an important feature of
  eukaryotic genomes and can bring regulatory sequences, such as enhancers or 
  transcription factor binding sites, in close physical proximity of regulated 
  target genes.
  Here, we provide 7C, a tool that uses protein binding signals from ChIP-seq and
  sequence motif information to predict chromatin looping events. Cross-linking 
  of proteins that bind close to loop anchors result in ChIP-seq signals at both 
  anchor loci. These signals are used at CTCF motif pairs together with their 
  distance and orientation to each other to predict whether they interact or 
  not. The resulting chromatin loops can be used to associate enhancers or 
  transcription factor binding sites (e.g. ChIP-seq peaks) to regulated target 
  genes.
  
## Intallation

The *sevenC* package depends on some R packages from Bioconductor. They can be 
installed as follows

```R
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer", "InteractionSet")
```

The development version can be installed from github using 
devtools:

```R
#install.packages("devtools")
devtools::install_github("ibn-salem/sevenC")
```

## Basic usage example
Here we show how to use the package to predict chromatin looping interactions 
among CTCF moif locations on chromosome 22. 
As input only a single bigWig file is used from a STAT1 ChIP-seq experiment 
in human GM12878 cells. 

#### Get motif pairs
```R
library(sevenC)

# load provided CTCF motifs in human genome
motifs <- motif.hg19.CTCF.chr22

# get motifs pairs
gi <- prepareCisPairs(motifs, maxDist = 10^6, scoreColname = "sig")
```

#### Add ChIP-seq data and compute correaltion
```R

# use example ChIP-seq bigWig file
bigWigFile <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig", 
  package = "sevenC")

# add ChIP-seq coverage and compute correaltion at motif pairs
gi <- addCor(gi, bigWigFile)
```

####  Predict loops

```R
# predict looping interactions among all motif pairs
loops <- predLoops(gi)
```

For more detailed usage instructions, see the package 
[vignette](https://ibn-salem.github.io/sevenC/articles/sevenC.html) or 
[reference documentation](https://ibn-salem.github.io/sevenC/reference/index.html).


## Issues
Please report any issues here: https://github.com/ibn-salem/sevenC/issues
