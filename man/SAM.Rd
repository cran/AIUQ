% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAM.R
\name{SAM}
\alias{SAM}
\title{Scattering analysis of microscopy}
\usage{
SAM(
  intensity = NA,
  intensity_str = "T_SS_mat",
  pxsz = 1,
  sz = c(NA, NA),
  mindt = 1,
  AIUQ_thr = c(1, 1),
  model_name = "BM",
  sigma_0_2_ini = NaN,
  param_initial = NA,
  num_optim = 1,
  msd_fn = NA,
  msd_grad_fn = NA,
  num_param = NA,
  uncertainty = FALSE,
  M = 50,
  sim_object = NA,
  msd_truth = NA,
  method = "AIUQ",
  index_q_AIUQ = NA,
  index_q_DDM = NA,
  message_out = TRUE,
  A_neg = "abs",
  square = FALSE,
  output_dqt = FALSE,
  output_isf = FALSE,
  output_modeled_isf = FALSE,
  output_modeled_dqt = FALSE
)
}
\arguments{
\item{intensity}{intensity profile. See 'Details'.}

\item{intensity_str}{structure of the intensity profile, options from
('SST_array','S_ST_mat','T_SS_mat'). See 'Details'.}

\item{pxsz}{size of one pixel in unit of micron, 1 for simulated data}

\item{sz}{frame size of the intensity profile in x and y directions,
number of pixels contained in each frame equals sz_x by sz_y.}

\item{mindt}{minimum lag time, 1 for simulated data}

\item{AIUQ_thr}{threshold for wave number selection, numeric vector of two
elements with values between 0 and 1. See 'Details'.}

\item{model_name}{fitted model, options from  ('BM','OU','FBM','OU+FBM',
'user_defined'), with Brownian motion as the default model. See 'Details'.}

\item{sigma_0_2_ini}{initial value for background noise. If NA, use minimum
value of absolute square of intensity profile in reciprocal space.}

\item{param_initial}{initial values for param estimation.}

\item{num_optim}{number of optimization.}

\item{msd_fn}{user defined mean squared displacement(MSD) structure, a
function of parameters and lag times. NA if \code{model_name} is not
'user_defined'.}

\item{msd_grad_fn}{gradient for user defined mean squared displacement
structure. If \code{NA}, then numerical gradient will be used for parameter
estimation in \code{'user_defined'} model.}

\item{num_param}{number of parameters need to be estimated in the intermediate
scattering function, need to be non-NA value for user_defined' model.}

\item{uncertainty}{a logical evaluating to TRUE or FALSE indicating whether
parameter uncertainty should be computed.}

\item{M}{number of particles. See 'Details'.}

\item{sim_object}{NA or an S4 object of class \code{simulation}.}

\item{msd_truth}{true MSD or reference MSD value.}

\item{method}{methods for parameter estimation, options from ('AIUQ','DDM_fixedAB','DDM_estAB').}

\item{index_q_AIUQ}{index range for wave number when using AIUQ method. See 'Details'.}

\item{index_q_DDM}{index range for wave number when using DDM method. See 'Details'.}

\item{message_out}{a logical evaluating to TRUE or FALSE indicating whether
or not to output the message.}

\item{A_neg}{controls modification for negative A(q), options from  ('abs','zero'),
with setting negative A(q) to its absolute value as the default.}

\item{square}{a logical evaluating to TRUE or FALSE indicating whether or not
to crop the original intensity profile into square image.}

\item{output_dqt}{a logical evaluating to TRUE or FALSE indicating whether or
not to compute observed dynamic image structure function(Dqt).}

\item{output_isf}{a logical evaluating to TRUE or FALSE indicating whether or
not to compute empirical intermediate scattering function(ISF).}

\item{output_modeled_isf}{a logical evaluating to TRUE or FALSE indicating
whether or not to compute modeled intermediate scattering function(ISF).}

\item{output_modeled_dqt}{a logical evaluating to TRUE or FALSE indicating
whether or not to compute modeled dynamic image structure function(Dqt).}
}
\value{
Returns an S4 object of class \code{SAM}.
}
\description{
Fast parameter estimation in scattering analysis of microscopy, using either
AIUQ or DDM method.
}
\details{
For simulated data using \code{simulation} in AIUQ package, \code{intensity}
will be automatically extracted from \code{simulation} class.

By default \code{intensity_str} is set to 'T_SS_mat', a time by space\eqn{\times}{\%\times}space
matrix, which is the structure of intensity profile obtained from \code{simulation}
class. For \code{intensity_str='SST_array'} , input intensity profile should be a
space by space by time array, which is the structure from loading a tif file.
For \code{intensity_str='S_ST_mat'}, input intensity profile should be a
space by space\eqn{\times}{\%\times}time matrix.

By default \code{AIUQ_thr} is set to \code{c(1,1)}, uses information from all
complete q rings. The first element affects maximum wave number selected,
and second element controls minimum proportion of wave number selected. By
setting 1 for the second element, if maximum wave number selected is less
than the wave number length, then maximum wave number selected is coerced to
use all wave number unless user defined another index range through \code{index_q_AIUQ}.

If \code{model_name} equals 'user_defined', or NA (will coerced to
'user_defined'), then \code{msd_fn} and \code{num_param} need to be provided
for parameter estimation.

Number of particles \code{M} is set to 50 or automatically extracted from
\code{simulation} class for simulated data using \code{simulation} in AIUQ
package.

By default, using all wave vectors from complete q ring, unless user defined
index range through \code{index_q_AIUQ} or \code{index_q_DDM}.
}
\examples{
library(AIUQ)
# Example 1: Estimation for simulated data
sim_bm = simulation(len_t=100,sz=100,sigma_bm=0.5)
show(sim_bm)
sam = SAM(sim_object = sim_bm)
show(sam)
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
