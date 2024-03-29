% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{strg2vcf}
\alias{strg2vcf}
\title{Create vcf string for a single SNP}
\usage{
strg2vcf(reference, alternative, total)
}
\arguments{
\item{reference}{an integer representing the number of reads with the
reference allele}

\item{alternative}{an integer representing the number of reads with the
alternative allele}

\item{total}{an integer representing the total number of reads observed at
this SNP.}
}
\value{
a character string coded as R,A:DP.
}
\description{
Creates a string with the information for a single SNP. The information is
coded as R,A:DP. R is the number of reads of the reference allele, A is the
number of reads of the alternative allele and DP is the total depth of
coverage.
}
\keyword{internal}
