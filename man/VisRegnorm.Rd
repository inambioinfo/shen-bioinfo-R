\name{VisRegnorm}
\alias{VisRegnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize gene count normalization based on regression.
%%  ~~function to do ... ~~
}
\description{
    Loess and size factor normalization from DESeq package are used. This 
    function is used to visualze the fit and do comparisons.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VisRegnorm(gene, spikein, x.col=1, y.col=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{Gene count data as a matrix or data.frame.}
  \item{spikein}{Spikein count data as a matrix or data.frame.}
  \item{x.col}{X sample column position or sample name. Can be integer or 
  character.}
  \item{y.col}{Y sample column position or sample name. Can be integer or 
  character.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
