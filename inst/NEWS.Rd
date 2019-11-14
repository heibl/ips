\name{ips-news}
\alias{ips-news} 
\title{News and changes to the ips (formerly: phyloch) package}
\section{VERSION 0.0.12}{
\subsection{NEW FEATURES}{
  \itemize{
  \item The new function \code{blastn} provides an interface to the command-line program blastn for nucleotide to nucleotide comparison.
  }
   \item The new function \code{combMyTree} grafts polytomies on a phylogeny.
  }
}
}
\section{VERSION 0.0.11}{
\subsection{NEW FEATURES}{
  \itemize{
  \item The new function \code{forceEqualTipHeights} corrects small rounding errors in edge lengths such that the resulting phylogeny will pass \code{ape::is.ultrametric}.
  }
}
\subsection{IMPROVED FEATURES}{
  \itemize{
  \item The new argument \code{weights} in \code{raxml} allows to assign individual weights to each column of the alignment. It corresponds to the -a flag in RAxML (see The RAxML v8.2.x Manual for details)
  }
}
}
\section{VERSION 0.0-10}{
\subsection{IMPROVED FEATURES}{
  \itemize{
  \item The new argument \code{exec} in \code{mrbayes} allows to specify the name and path of the MrBayes executable explicitly. If the executable is in the search path, \code{exec} can be missing.
  }
}
}
\section{VERSION 0.0-9}{
\subsection{IMPROVED FEATURES}{
  \itemize{
  \item \code{write.nex} (and \code{matrixBlocks}) were extended to handle standard (morphological, etc.) data in a data frame. This feature including the coding of ambiguous characters was tested successfully with MrBayes.
  }
}
\subsection{DEPRECATED FEATURES}{
   \itemize{
    \item \code{write.partioned.nex} was removed from the package; its functionality has been integrated into \code{write.nex}.
   } 
}
}
\section{VERSION 0.0-8}{
\subsection{IMPROVED FEATURES}{
   \itemize{
    \item \code{write.nex} can now handle multiple DNA sequence alignments.
    \item Argument \code{interleave} of functions \code{write.fas}, \code{write.phy}, and \code{write.nex} has been renamed to \code{block.width} for clarity.
   } 
  }
  \subsection{DEPRECATED FEATURES}{
   \itemize{
    \item \code{write.partioned.nex} will be removed soon from the package; its functionality has been integrated into \code{write.nex}.
   } 
  }
  \subsection{BUG FIXES}{
   \itemize{
    \item The argument \code{run = TRUE} in functions \code{mrbayes} and \code{mrbayes.mixed} was broken on Windows platforms. (Thanks to Liam Revell and Klaus Schliep for report ans fix).
   } 
  }
}
\section{VERSION 0.0-7}{ 
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{mafft} received the additional argument \code{options}, which can be used to request options such as e.g. --adjustdirection that are not build into the function's interface.
      \item This version includes a new internal function \code{phylo2mafft}, which does exactly the same thing as the RUBY script newick2mafft.rb on the MAFFT website (\url{http://mafft.cbrc.jp/alignment/software/newick2mafft.rb}): it converts a user-defined guide tree into a format readible by MAFFT.
      \item The functions \code{delete.empty.cells} and \code{fillEndsWithN} are now using Emanuel Paradis' bit-level coding for DNA sequences, which makes them much faster.
    }
  }
  \subsection{DEPRECATED FEATURES}{
    \itemize{
      \item The function \code{c.genes} has been superseeded by the cbind method for object of class \code{"DNAbin"} provided in the package ape; \code{c.genes} will be removed in one of the following versions.
    }
  }
}
