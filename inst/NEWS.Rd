\name{NEWS}
\title{Rchoice News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.3-6}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Adding the Hessian to ivplm.
      \item Fixing bug with the LR test in "hetprob()"
    }
  }  
}

\section{CHANGES IN VERSION 0.3-5}{
\subsection{NEW FEATURES}{
    \itemize{
      \item New urls in Description. 
      \item New README.
      \item New methods for effect.
      \item Improved documentation.
      \item AIC and BIC for different classes are dropped. We just need "logLik()".
      \item Changing the behavior of set.seed as a good practice. 
    }
  }  
}

\section{CHANGES IN VERSION 0.3-4}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Fixing bug with mtable.
    }
  }  
}

\section{CHANGES IN VERSION 0.3-3}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Add hetprob function to estimate heterokedastic binary models.
      \item Add ivplm function to estimatate instrumental variables probit model by Maximum Likelihood.
      \item Fix minor bug in the \code{model.frame.rFormula} function with new update of \pkg{plm}. 
      \item Slightly improved documentation. 
    }
  }  
}

\section{CHANGES IN VERSION 0.3-2}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Fixing errors and warnings
    }
  }  
}

\section{CHANGES IN VERSION 0.3-1}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Including CITATION
    }
  }  
}

\section{CHANGES IN VERSION 0.3}{
  \subsection{FIXED BUGS}{
    \itemize{
      \item fixed minor bug in the \code{model.frame.rFormula} function when the id variable is specified using \code{pdata.frame} function
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item the package not longer includes a vignette. The documentation can be dowloaded at \url{http://msarrias.weebly.com/rchoice-package-in-r.html}
      \item fix some notes regarding to the AIC and BIC methods
      \item statistic in \code{summary} is now labeled as z-value instead of t-value
      \item \code{cov.Rchoice}, \code{cor.Rchoice} and \code{se.cov.Rchoice} functions are now deprecated. All these functions are now available in the generic function \code{vcov}
      \item slightly improved documentation 
    }
  }
}

\section{CHANGES IN VERSION 0.2}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item  panel data or longitudinal data is now allow
      \item  added the function \code{se.cov.Rchoice} to compute the standard error of the variance-covariance matrix of random parameters
      \item more flexibility to model hirarchical variables. The argument \code{mvar} can be used to choose which variables enter in the mean of the random parameters.
      \item added the Johnson's Sb function for random parameters.
      \item added the \code{rFormula} and related functions
      \item initial values for standard deviation of random parameters are set to 0.1. Previously was 0. The argument \code{init.ran} can be used to set the initial values. For example \code{init.ran = 0} will set all the initial values for the standard deviation paramters at 0
      \item slightly improved documentation
    }
  }
}
