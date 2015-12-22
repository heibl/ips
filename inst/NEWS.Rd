\name{ips-news}
\alias{ips-news} 
\title{News and changes to the ips (formerly: phyloch) package}
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
\section{CHANGES IN PHYLOCH VERSION 1.5-6}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item There is a new function \code{fixedFontBinom}, inspired by Emmanuel Paradis' \code{fixedFontLabel} function in ape. The function turns the tiplabels of a tree into expressions and so can accomodate for taxonomic conventions, e.g., italic font for genera and species, the 'times' symbol for hybrid taxa and so on.
    }
  }
}
\section{CHANGES IN PHYLOCH VERSION 1.5-5}{
\itemize{
      \item Minor changes and recompilation for R 3.0.1.
    }
}
\section{CHANGES IN PHYLOCH VERSION 1.5-4}{
\itemize{
      \item Minor changes.
    }
}
\section{CHANGES IN PHYLOCH VERSION 1.5-3}{
\itemize{
      \item Minor changes and recompilation for R 3.0.0.
    }
}
\section{CHANGES IN PHYLOCH VERSION 1.5-2}{
  \subsection{IMPROVED FEATURES}{
  \itemize{
      \item \code{box.clades} and \code{box.tips} have been rewritten and work more cleanly now. Yet, they still support only directions  'rightwards' and 'leftwards' in phylograms. The argument 'align' in box.tips() has now an additional possible value ("tipwise"), which now is the default.
    }
  }
}
\section{CHANGES IN PHYLOCH VERSION 1.5-1}{
  \subsection{NEW FEATURES}{
  \itemize{
      \item Argument \code{quiet} added to \code{mafft} and \code{c.genes} to be able to suppress screen output by both functions.
    }
  }
}   
\section{CHANGES IN PHYLOCH VERSION 1.5-0}{
    \subsection{BUG FIXES}{
    \itemize{
      \item \code{fixNodes} did not work on some ladderized trees. This is (hopefully) now fixed, but users are still encouraged to report any peculiarities. A this stage it might still be good to doublecheck if a unknown tree is manipulated for the first time. (Thanks to Ingo Michalak and Johan Nylander for bug report.) 

    }
    }
  \subsection{NEW FEATURES}{
    \itemize{
      \item This ChangeLog is included.
      \item The stratigraphic chart 2012 is provided for \code{axisGeo}; see ?strat2012 for details. (Thanks to Lars Nauheimer for providing the data frame).
    }
  }
}