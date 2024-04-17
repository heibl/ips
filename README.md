
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ips: interfaces to phylogenetic software

> NOTE: `ips` contains most functions that were formerly included in the
> `phyloch` package plus some more recent functions. Those `phyloch`
> functions related to tree plotting have been moved to the
> [`viper`](https://github.com/heibl/viper) package.

This is a bundle of functions that present **i**nterfaces to popular
**p**hylogenetic **s**oftware for sequence alignment, masking of
sequence alignments, and estimation of phylogenies and ancestral
character states. In additions, there functions for reading,
manipulating and writing phylogenetic data (multiple sequence alignments
and phylogenetic trees).

## Introduction

There are several functions for reading and writing DNA sequences in
FASTA, PHYLIP, and NEXUS format: `read.fas`, `read.phy`, `read.nex`,
`write.fas`, `write.phy`, and `write.nex`. Some functions are available
for integrating BEAST with R. XML input files for BEAST can be generated
with `rbeauti`. Two functions are designed to read TreeAnnotator output:
`read.beast` will render an object of class `phylo` with additional node
statistics appended as list elements. These additional node statistics
will be lost be the subsequent use of `ladderize` or `rotate` (or
similar functions that change the ordering of internal
nodes).`read.beast.table` also parses the TreeAnnotator output, but
returns a matrix of node statistics. This package itself does not
implement techniques for phylogenetic analyses, but provides a series of
wrappers for commonly used software packages. Sequence alignment can be
done with the `mafft` and `prank`; cleaning of sequences with `gblocks`
and `aliscore`. The function `raxml` and `mrbayes` are intended for
phylogenetic tree search. Running `mrbayes` with argument `run = FALSE`
can be used to create MrBayes-executable NEXUS files. Finally, wrappers
are provided for `Multistate` in the `BayesTraits` package (see
`multistateML` and `multistateMCMC`). Several plotting functions
(`HPDbars`, `clade.bars`, `box.clades`, `box.tips`, `tip.color`,
`edge.color` have been moved to the `viper` package.

## Installation

`ips` is available via CRAN and can be installed:

``` r
install.packages("ips")
```

Or you can install the development version via GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("heibl/ips")
```
