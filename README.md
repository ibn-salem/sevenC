# chromloop
Prediction of Chromatin Looping Interaction from ChIP-seq and Sequnece Motifs

## Intallation

```R
#install.packages("devtools")
devtools::install_github("ibn-salem/chromloop")
```

## Example usage

```R
# load CTCF motifs
motifs <- motif.hg19.CTCF

# add ChIP-seq coverage
motifs <- addCovToGR(motifs, "coverage.bigWig")

# get pairs of motifs
gi <- getCisPairs(motifs, maxDist = 10^6)

# compute correlation of profiles
gi <- applyToCloseGI(gi, "cov", fun = cor)

# add known loops
gi <- addInteractionSupport(gi, knownLoops)

# train model 
fit <- glm(loop ~ cov + dist + strandOrientation, data = mcols(gi), familiy = "logit")

# predict loops
predicted_loops <- predict(fit, type = "response")

```

