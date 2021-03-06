% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gof_lighthergm.R
\name{swap_formula_network}
\alias{swap_formula_network}
\title{Swaps the network on the lhs of a formula for a new one with the given environment}
\usage{
swap_formula_network(new_net, net_formula, env)
}
\arguments{
\item{new_net}{A network object to be inserted into the lhs of the formula}

\item{net_formula}{The target formula}

\item{env}{The environment to assign to the formula}
}
\value{
A new formula with the lhs swapped
}
\description{
Swaps the network on the lhs of a formula for a new one with the given environment
}
