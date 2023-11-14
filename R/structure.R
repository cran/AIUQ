#' SAM class
#'
#'@description
#' S4 class for fast parameter estimation in scattering analysis of microscopy,
#' using either \code{AIUQ} or \code{DDM} method.
#'
#' @slot pxsz numeric.  Size of one pixel in unit of micron with default value 1.
#' @slot mindt numeric. Minimum lag time with default value 1.
#' @slot sz numeric. Frame size of the intensity profile, number of pixels
#' contained in each frame equals sz^2.
#' @slot len_t integer. Number of time steps.
#' @slot len_q integer. Number of wave vector.
#' @slot q vector. Wave vector in unit of um^-1.
#' @slot d_input vector. Sequence of lag times.
#' @slot anisotropic logical. A logical evaluating to \code{TRUE} or \code{FALSE}
#' indicating whether anisotropic estimation should be performed
#' @slot B_est_ini numeric. Estimation of B. This parameter is determined by the
#' noise in the system. See 'References'.
#' @slot A_est_ini vector. Estimation of A(q). Note this parameter is
#' determined by the properties of the imaged material and imaging optics.
#' See 'References'.
#' @slot I_o_q_2_ori vector. Absolute square of Fourier transformed intensity
#' profile, ensemble over time.
#' @slot model_name character. Fitted model, options from
#' ('BM','OU','FBM','OU+FBM', 'user_defined').
#' @slot param_est vector. Estimated parameters contained in MSD.
#' @slot sigma_2_0_est numeric. Estimated variance of background noise.
#' @slot msd_est vector. Estimated MSD.
#' @slot uncertainty logical. A logical evaluating to TRUE or FALSE indicating whether
#' parameter uncertainty should be computed.
#' @slot msd_lower vector. Lower bound of 95% confidence interval of MSD.
#' @slot msd_upper vector. Upper bound of 95% confidence interval of MSD.
#' @slot msd_truth vector. True MSD or reference MSD value.
#' @slot sigma_2_0_truth vector.  True variance of background noise, non NA for
#' simulated data using \code{simulation}.
#' @slot param_truth vector. True parameters used to construct MSD, non NA for
#' simulated data using \code{simulation}.
#' @slot index_q vector. Selected index of wave vector.
#' @slot Dqt matrix. Dynamic image structure function D(q,delta t).
#' @slot I_q matrix. Fourier transformed intensity profile with structure 'SS_T_mat'.
#' @slot AIC numeric. Akaike information criterion score.
#' @slot mle numeric. Maximum log likelihood value.
#' @slot param_uq_range  matrix. 95% confidence interval for estimated parameters.
#'
#' @method show SAM
#' @author \packageAuthor{AIUQ}
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
#'
#' Cerbino, R., & Trappe, V. (2008). Differential dynamic microscopy: probing
#' wave vector dependent dynamics with a microscope. Physical review letters,
#' 100(18), 188102.
#'
#' @keywords classes
methods::setClass("SAM", representation(
  #intensity_str = "character",
  pxsz = "numeric",
  mindt = "numeric",
  sz = "numeric",
  len_t = "integer",
  len_q = "integer",
  q = "vector",
  d_input = "vector",
  anisotropic = "logical",
  B_est_ini = "numeric",
  A_est_ini = "vector",
  #num_q_max = "numeric",
  I_o_q_2_ori = "vector",
  #q_ori_ring_loc_unique_index = "list",
  model_name = "character",
  param_est = "vector",
  sigma_2_0_est = "numeric",
  msd_est = "vector",
  uncertainty  = "logical",
  msd_lower = "vector",
  msd_upper = "vector",
  msd_truth = "vector",
  sigma_2_0_truth = "vector",
  param_truth = "vector",
  method = "character",
  index_q = "vector",
  Dqt = "matrix",
  I_q = "matrix",
  AIC = "numeric",
  mle = "numeric",
  param_uq_range = "matrix"
  #p = "numeric"
)
)

## Show
if(!isGeneric("show")){
  setGeneric(name = "show",
             def = function(object) standardGeneric("show"))
}

setMethod("show", "SAM",
          function(object){show.sam(object)})

#' Simulation class
#'
#' @description
#' S4 class for 2D particle movement simulation.
#'
#' @slot sz integer. Frame size of the intensity profile, number of pixels
#' contained in each frame equals sz^2.
#' @slot len_t integer. Number of time steps.
#' @slot noise character. Background noise, options from ('uniform','gaussian').
#' @slot model_name character. Simulated stochastic process, options from ('BM','OU','FBM','OU+FBM').
#' @slot M integer. Number of particles.
#' @slot pxsz numeric.  Size of one pixel in unit of micron, 1 for simulated data.
#' @slot mindt numeric. Minimum lag time, 1 for simulated data.
#' @slot pos matrix. Position matrix for particle trajectory, see 'Details'.
#' @slot intensity matrix. Filled intensity profile, see 'Details'.
#' @slot num_msd vector. Numerical mean squared displacement (MSD).
#' @slot param vector. Parameters used to construct MSD.
#' @slot theor_msd vector. Theoretical MSD.
#' @slot sigma_2_0 vector. Variance of background noise.
#'
#' @method show simulation
#' @details
#' \code{intensity} should has structure 'T_SS_mat', matrix with dimension
#' \code{len_t} by \code{sz}\eqn{\times}{%\times}\code{sz}.
#'
#' \code{pos} should be the position matrix with dimension
#' \code{M}\eqn{\times}{%\times}\code{len_t}. See \code{\link{bm_particle_intensity}},
#' \code{\link{ou_particle_intensity}}, \code{\link{fbm_particle_intensity}},
#' \code{\link{fbm_ou_particle_intensity}}.
#'
#' @author \packageAuthor{AIUQ}
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
#'
#' Cerbino, R., & Trappe, V. (2008). Differential dynamic microscopy: probing
#' wave vector dependent dynamics with a microscope. Physical review letters,
#' 100(18), 188102.
#'
#'
#'@keywords classes
methods::setClass("simulation", representation(
  sz = "integer",
  len_t = "integer",
  noise = "character",
  model_name = "character",
  M = "integer",
  pxsz = "numeric",
  mindt = "numeric",
  pos = "matrix", #first M equals pos0
  intensity = "matrix",
  num_msd = "vector",
  param = "vector",
  theor_msd = "vector",
  sigma_2_0 = "vector"
)
)

## Show and plot
if(!isGeneric("show")){
  setGeneric(name = "show",
             def = function(object) standardGeneric("show"))
}

setMethod("show", "simulation",
          function(object){show.simulation(object)})






