\name{svynls}
\alias{svynls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability-weighted nonlinear least squares
}
\description{
Fits a nonlinear model by probability-weighted least squares.  Uses
\code{nls} to do the fitting, but estimates design-based standard errors with either
linearisation or replicate weights. See \code{\link{nls}} for 
documentation of model specification and fitting.
}
\usage{
svynls(formula, design, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    Nonlinear model specified as a formula; see \code{\link{nls}}
}
  \item{design}{
    Survey design object
  }
  \item{\dots}{
    Other arguments to \code{nls} (especially, \code{start}). Also
    supports \code{return.replicates} for replicate-weight designs and
    \code{influence} for other designs. 
  }
}

\value{
  Object of class \code{svynls}. The fitted \code{nls} object is
  included as the \code{fit} element.
}

\seealso{
\code{\link{svymle}} for maximum likelihood with linear predictors on
one or more parameters
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
