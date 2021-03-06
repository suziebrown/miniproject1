% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancestry_size.R
\name{ancestrySize}
\alias{ancestrySize}
\title{Sample ancestry size}
\usage{
ancestrySize(history, sampl = NULL, maxgen = NULL)
}
\arguments{
\item{history}{an ocject of class genealogy, as returned from WFsim or moranSim}

\item{sampl}{a vector of integers indicating which individuals from the most recent generation are to be included in the sample}

\item{maxgen}{the maximum number of generations back in time to investigate}
}
\value{
a vector indicating the number of ancestors in each generation 1,...,N.gen
}
\description{
Find the number of ancestors of a given sample in each generation
}
\author{
Suzie Brown
}
