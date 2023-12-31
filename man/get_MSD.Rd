% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get_MSD}
\alias{get_MSD}
\title{Construct MSD}
\usage{
get_MSD(theta, d_input, model_name, msd_fn = NA)
}
\arguments{
\item{theta}{parameters in MSD function}

\item{d_input}{sequence of lag times}

\item{model_name}{model name for the process, options from ('BM','OU','FBM',
'OU+FBM','user_defined')}

\item{msd_fn}{user defined mean squared displacement structure (MSD), a
function of \code{param} parameters and \code{d_input} lag times}
}
\value{
A vector of MSD values for a given sequence of lag times.
}
\description{
Construct estimated mean squared displacement (MSD) for a given stochastic process.
}
\details{
For Brownian Motion, the MSD follows
\deqn{MSD_{BM}(\Delta t) = \theta_1\Delta t= 4D\Delta t}{\%MSD_{BM}(\Delta t) = \theta_1\Delta t= 4D\Delta t}
where \code{D} is the diffusion coefficient.

For Ornstein–Uhlenbeck process,  the MSD follows
\deqn{MSD_{OU}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})}{\%MSD_{OU}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})}
where \eqn{\theta_1=\rho}{\%\theta_1=\rho}
is the correlation with previous steps.

For fractional Brownian Motion,  the MSD follows
\deqn{MSD_{FBM}(\Delta t) =\theta_1\Delta t^{\theta_2}}{\%MSD_{FBM}(\Delta t) =\theta_1\Delta t^{\theta_2}}
where \eqn{\theta_2=2H}{\%\theta_2=2H} with \code{H} is the the Hurst parameter.

For 'OU+FBM', the MSD follows
\deqn{MSD_{OU+FBM}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})+\theta_3\Delta t^{\theta_4}}{\%MSD_{OU+FBM}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})+\theta_3\Delta t^{\theta_4}}
}
\examples{
library(AIUQ)
# Construct MSD for BM
get_MSD(theta=0.2,d_input=0:100,model_name='BM')

}
\references{
Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
quantification in scattering analysis of microscopy.
arXiv preprint arXiv:2309.02468.

Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
Uncertainty quantification and estimation in differential dynamic microscopy.
Physical Review E, 104(3), 034610.
}
\author{
\packageAuthor{AIUQ}
}
\keyword{internal}
