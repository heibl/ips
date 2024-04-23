# ips 0.0.12

## New features

* `blastn()` provides an interface to the command-line program blastn for nucleotide to nucleotide comparison.
  
* `combMyTree()` grafts polytomies on a phylogeny.
   
* `deleteEmptyCells()` has now a sister function `identifyEmtpyCells()` that won't delete but only identify the empty rows and columns in a sequence alignment. This is useful if there is a second matrix (e.g. with confidence scores) tied to the sequence alignment.

# ips 0.0.11

## New features

* `forceEqualTipHeights()` corrects small rounding errors in edge lengths such that the resulting phylogeny will pass `ape::is.ultrametric()`.

## Improved features

* Argument `weights` in `raxml()` allows to assign individual weights to each column of the alignment. It corresponds to the -a flag in RAxML (see The RAxML v8.2.x Manual for details)

# ips 0.0-10

## Improved features

* Argument `exec` in `mrbayes()` allows to specify the name and path of the MrBayes executable explicitly. If the executable is in the search path, `exec` can be missing.

# ips 0.0-9

## Improved features

* `write.nex()` and `matrixBlock` were extended to handle standard (morphological, etc.) data in a data frame. This feature including the coding of ambiguous characters was tested successfully with MrBayes.

## Deprecated features

* `write.partioned.nex()` was removed from the package; its functionality has been integrated into `write.nex()`.

# ips 0.0-8

## Improved features

* `write.nex()` can now handle multiple DNA sequence alignments.

* Argument `interleave` of functions `write.fas()`, `write.phy()`, and `write.nex()` has been renamed to `block.width` for clarity.

## Deprecated features

* `write.partioned.nex()` will be removed soon from the package; its functionality has been integrated into `write.nex()`.

## Bug fixes

 * Argument `run = TRUE` in `mrbayes()` and `mrbayes.mixed()` was broken on Windows platforms. (Thanks to Liam Revell and Klaus Schliep for report and fix).

# ips 0.0-7

## New features

* `mafft()` received the additional argument `options`, which can be used to request options such as e.g.  `--adjustdirection` that are not build into the function's interface.

* This version includes a new internal function `phylo2mafft()`, which does exactly the same thing as the RUBY script `newick2mafft.rb` on the MAFFT website (\url{https://mafft.cbrc.jp/alignment/software/newick2mafft.rb}): it converts a user-defined guide tree into a format readible by MAFFT.

* `delete.empty.cells()` and `fillEndsWithN()` are now using Emanuel Paradis' bit-level coding for DNA sequences, which makes them much faster.

## Deprecated features

* `c.genes()` has been superseeded by the cbind method for objects of class `"DNAbin"` provided in the package **ape**; `c.genes()` will be removed in one of the following versions.
