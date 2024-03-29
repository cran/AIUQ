% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{SAM_Dqt}
\alias{SAM_Dqt}
\title{Compute dynamic image structure function}
\usage{
SAM_Dqt(len_q, index_q, len_t, I_q_matrix, q_ori_ring_loc_unique_index, sz)
}
\arguments{
\item{len_q}{number of wave number}

\item{index_q}{a vector of selected wave number index}

\item{len_t}{number of time steps}

\item{I_q_matrix}{intensity profile in reciprocal space (after Fourier transformation)}

\item{q_ori_ring_loc_unique_index}{index for wave vector that give unique frequency}

\item{sz}{frame size of intensity profile}
}
\value{
Matrix of dynamic image structure with dimension \code{len_q} by \code{len_t-1}.
}
\description{
Compute dynamic image structure function(Dqt) using Fourier transformed
intensity profile and a selection of wave number(q) range.
}
\details{
Dynamic image structure function(Dqt) can be obtained from ensemble average
of absolute values squared of Four transformed intensity difference:
\deqn{D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}{\%D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}
See 'References'.
}
\references{
Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
quantification in scattering analysis of microscopy.
arXiv preprint arXiv:2309.02468.

Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
Uncertainty quantification and estimation in differential dynamic microscopy.
Physical Review E, 104(3), 034610.

Cerbino, R., & Trappe, V. (2008). Differential dynamic microscopy: probing
wave vector dependent dynamics with a microscope. Physical review letters,
100(18), 188102.
}
\author{
\packageAuthor{AIUQ}
}
\keyword{internal}
