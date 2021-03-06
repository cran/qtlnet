\name{loci.qtlnet}
\alias{loci.qtlnet}
\alias{est.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QTL architecture per node as list}
\description{
Determines QTL that affect each phenotype conditional on the model-averaged
network and on covariates.
}
\usage{
loci.qtlnet(qtlnet.object, chr.pos = TRUE, merge.qtl = 10, \dots)
est.qtlnet(qtlnet.object, ..., verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qtlnet.object}{Object of class \code{qtlnet}.}
  \item{chr.pos}{Include chromsome and position if \code{TRUE}.}
  \item{merge.qtl}{Merge QTL within \code{merge.qtl} cM of the mean QTL
    per chromosome across all nodes. No merge if 0 or less. This can
    reduce the number of QTL nodes to one per chr.}
  \item{\dots}{Additional unused arguments.}
  \item{verbose}{verbose output if \code{TRUE}.}
}
\value{
List containing, for each phenotype in the network, a character vector
of the QTL names as \code{chr@pos}, or pseudomarker name if
\code{chr.pos} is \code{FALSE}.
}
\author{
Brian S. Yandell and Elias Chaibub Neto
}
\references{
Chaibub Neto E, Keller MP, Attie AD, Yandell BS (2010)
Causal Graphical Models in Systems Genetics: a unified
framework for joint inference of causal network and
genetic architecture for correlated phenotypes.
Ann Appl Statist 4: 320-339.
\url{http://www.stat.wisc.edu/~yandell/doc/2010/92.AnnApplStat.pdf}
}
\seealso{
\code{\link{mcmc.qtlnet}}
}
\examples{
loci.qtlnet(Pscdbp.qtlnet)
}

