## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=6, fig.height=4)
#, global.par=TRUE
options(digits = 4)

## -----------------------------------------------------------------------------
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
tree <- pratchet(primates, trace=0) |> acctran(primates)
parsimony(tree, primates)

## -----------------------------------------------------------------------------
anc.acctran <- ancestral.pars(tree, primates, "ACCTRAN")
anc.mpr <- ancestral.pars(tree, primates, "MPR")

## ----seqLogo, fig.cap="Fig 1. Ancestral reconstruction for a node.", eval=FALSE----
#  library(seqLogo)
#  seqLogo( t(subset(anc.mpr, getRoot(tree), 1:20)[[1]]), ic.scale=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("seqLogo")

## ----MPR, fig.cap="Fig 2. Ancestral reconstruction using MPR."----------------
plotAnc(tree, anc.mpr, 17)
title("MPR")

## ----ACCTRAN, fig.cap="Fig 3. Ancestral reconstruction using ACCTRAN."--------
plotAnc(tree, anc.acctran, 17)
title("ACCTRAN")

## -----------------------------------------------------------------------------
fit <- pml(tree, primates)
fit <- optim.pml(fit, model="F81", control = pml.control(trace=0))

## -----------------------------------------------------------------------------
anc.ml <- ancestral.pml(fit, "ml")
anc.bayes <- ancestral.pml(fit, "bayes")

## ----plotML, fig.cap="Fig 4. Ancestral reconstruction the using the maximum likelihood."----
plotAnc(tree, anc.ml, 17)
title("ML")

## ----plotB, fig.cap="Fig 5. Ancestral reconstruction using (empirical) Bayes."----
plotAnc(tree, anc.bayes, 17)
title("Bayes")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

