# chromloop
Prediction of Chromatin Looping Interaction from ChIP-seq and Sequnece Motifs

## Intallation

```R
#install.packages("devtools")
devtools::install_github("ibn-salem/chromloop")
```

## Usage example
Here we show how to use the package to predict chromatin looping interactions 
among CTCF moif locations on chromosome 22. 
As input only a single coverage file is used from a STAT1 ChIP-seq experiment 
in human GM12878 cells. 

### Load motifs and ChIP-seq data
```R
library(chromloop)

# load provided CTCF motifs
motifs <- motif.hg19.CTCF.chr22

# use example ChIP-seq coverage file
bigWigFile <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig", 
  package = "chromloop")

# add ChIP-seq coverage
motifs <- addCovToGR(motifs, bigWigFile)

```

### Get pairs and correlation of ChIP-seq coverage
```R
# get pairs of motifs
gi <- getCisPairs(motifs, maxDist = 10^6)

# add motif orientation
gi <- addStrandCombination(gi)

# compute correlation of ChIP-seq profiles
gi <- applyToCloseGI(gi, "cov", fun = cor)
```

### Train predictor with known loops
```R
# parse known loops
knownLoopFile <- system.file("extdata", 
  "GM12878_HiCCUPS.chr22_1-18000000.loop.txt", package="chromloop")

knownLoops <- parseLoopsRao(knownLoopFile)

# add known loops
gi <- addInteractionSupport(gi, knownLoops)

# train model 
fit <- glm(loop ~ cor + dist + strandOrientation, 
  data = mcols(gi), 
  family = binomial("logit"))
```

### Predict loops
```R
# predict loops
gi$pred <- predict(fit, type = "response", newdata = mcols(gi)) 

# plot prediction score 
boxplot(gi$pred ~ gi$loop)

```

