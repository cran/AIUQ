#' SAM class
#'
#'@description
#' S4 class for fast parameter estimation in scattering analysis of microscopy,
#' using either \code{AIUQ} or \code{DDM} method.
#'
#' @slot pxsz numeric.  Size of one pixel in unit of micron with default value 1.
#' @slot mindt numeric. Minimum lag time with default value 1.
#' @slot sz vector. Frame size of the intensity profile in x and y directions,
#' number of pixels contained in each frame equals sz_x by sz_y.
#' @slot len_t integer. Number of time steps.
#' @slot len_q integer. Number of wave vector.
#' @slot q vector. Wave vector in unit of um^-1.
#' @slot d_input vector. Sequence of lag times.
#' @slot B_est_ini numeric. Estimation of B. This parameter is determined by the
#' noise in the system. See 'References'.
#' @slot A_est_ini vector. Estimation of A(q). Note this parameter is
#' determined by the properties of the imaged material and imaging optics.
#' See 'References'.
#' @slot I_o_q_2_ori vector. Absolute square of Fourier transformed intensity
#' profile, ensemble over time.
#' @slot q_ori_ring_loc_unique_index list. List of location index of non-duplicate
#' values for each q ring.
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
#' @slot ISF matrix. Empirical intermediate scattering function f(q,delta t).
#' @slot I_q matrix. Fourier transformed intensity profile with structure 'SS_T_mat'.
#' @slot AIC numeric. Akaike information criterion score.
#' @slot mle numeric. Maximum log likelihood value.
#' @slot param_uq_range  matrix. 95% confidence interval for estimated parameters.
#' @slot modeled_Dqt matrix. Modeled dynamic image structure function D(q,delta t).
#' @slot modeled_ISF matrix. Modeled intermediate scattering function f(q,delta t).
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
  sz = "vector",
  len_t = "integer",
  len_q = "integer",
  q = "vector",
  d_input = "vector",
  B_est_ini = "numeric",
  A_est_ini = "vector",
  #num_q_max = "numeric",
  I_o_q_2_ori = "vector",
  q_ori_ring_loc_unique_index = "list",
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
  ISF = "matrix",
  I_q = "matrix",
  AIC = "numeric",
  mle = "numeric",
  param_uq_range = "matrix",
  modeled_Dqt = "matrix",
  modeled_ISF = "matrix"
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

#' Anisotropic SAM class
#'
#'@description
#' S4 class for fast parameter estimation in scattering analysis of microscopy
#' for anisotropic processes, using either \code{AIUQ} or \code{DDM} method.
#'
#' @slot pxsz numeric.  Size of one pixel in unit of micron with default value 1.
#' @slot mindt numeric. Minimum lag time with default value 1.
#' @slot sz vector. Frame size of the intensity profile in x and y directions,
#' number of pixels contained in each frame equals sz_x by sz_y.
#' @slot len_t integer. Number of time steps.
#' @slot len_q integer. Number of wave vector.
#' @slot q vector. Wave vector in unit of um^-1.
#' @slot d_input vector. Sequence of lag times.
#' @slot B_est_ini numeric. Estimation of B. This parameter is determined by the
#' noise in the system. See 'References'.
#' @slot A_est_ini vector. Estimation of A(q). Note this parameter is
#' determined by the properties of the imaged material and imaging optics.
#' See 'References'.
#' @slot I_o_q_2_ori vector. Absolute square of Fourier transformed intensity
#' profile, ensemble over time.
#' @slot q_ori_ring_loc_unique_index list. List of location index of non-duplicate
#' values for each q ring.
#' @slot model_name character. Fitted model, options from
#' ('BM','OU','FBM','OU+FBM', 'user_defined').
#' @slot param_est matrix. Estimated parameters contained in MSD.
#' @slot sigma_2_0_est vector. Estimated variance of background noise.
#' @slot msd_est matrix. Estimated MSD.
#' @slot uncertainty logical. A logical evaluating to TRUE or FALSE indicating whether
#' parameter uncertainty should be computed.
#' @slot msd_truth matrix. True MSD or reference MSD value.
#' @slot sigma_2_0_truth vector.  True variance of background noise, non NA for
#' simulated data using \code{simulation}.
#' @slot param_truth matrix. True parameters used to construct MSD, non NA for
#' simulated data using \code{aniso_simulation}.
#' @slot index_q vector. Selected index of wave vector.
#' @slot I_q matrix. Fourier transformed intensity profile with structure 'SS_T_mat'.
#' @slot AIC numeric. Akaike information criterion score.
#' @slot mle numeric. Maximum log likelihood value.
#' @slot msd_x_lower vector. Lower bound of 95% confidence interval of MSD in x directions.
#' @slot msd_x_upper vector. Upper bound of 95% confidence interval of MSD in x directions.
#' @slot msd_y_lower vector. Lower bound of 95% confidence interval of MSD in y directions.
#' @slot msd_y_upper vector. Upper bound of 95% confidence interval of MSD in y directions.
#' @slot param_uq_range  matrix. 95% confidence interval for estimated parameters.
#'
#' @method show aniso_SAM
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
methods::setClass("aniso_SAM", representation(
  pxsz = "numeric",
  mindt = "numeric",
  sz = "vector",
  len_t = "integer",
  len_q = "integer",
  q = "vector",
  d_input = "vector",
  B_est_ini = "numeric",
  A_est_ini = "vector",
  #num_q_max = "numeric",
  I_o_q_2_ori = "vector",
  q_ori_ring_loc_unique_index = "list",
  model_name = "character",
  param_est = "matrix",
  sigma_2_0_est = "vector",
  msd_est = "matrix",
  uncertainty = "logical",
  msd_x_lower = "vector",
  msd_x_upper = "vector",
  msd_y_lower = "vector",
  msd_y_upper = "vector",
  msd_truth = "matrix",
  sigma_2_0_truth = "vector",
  param_truth = "matrix",
  method = "character",
  index_q = "vector",
  I_q = "matrix",
  AIC = "numeric",
  mle = "numeric",
  param_uq_range = "matrix"
)
)

## Show
if(!isGeneric("show")){
  setGeneric(name = "show",
             def = function(object) standardGeneric("show"))
}

setMethod("show", "aniso_SAM",
          function(object){show.aniso_sam(object)})

#' Simulation class
#'
#' @description
#' S4 class for 2D particle movement simulation.
#'
#' @slot sz vector. Frame size of the intensity profile, number of pixels
#' contained in each frame equals \code{sz[1]} by \code{sz[2]}.
#' @slot len_t integer. Number of time steps.
#' @slot noise character. Background noise, options from ('uniform','gaussian').
#' @slot model_name character. Simulated stochastic process, options from ('BM','OU','FBM','OU+FBM').
#' @slot M integer. Number of particles.
#' @slot pxsz numeric.  Size of one pixel in unit of micron, 1 for simulated data.
#' @slot mindt numeric. Minimum lag time, 1 for simulated data.
#' @slot pos matrix. Position matrix for particle trajectory, see 'Details'.
#' @slot intensity matrix. Filled intensity profile, see 'Details'.
#' @slot num_msd vector. Numerical mean squared displacement (MSD).
#' @slot param vector. Parameters for simulated stochastic process.
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
  sz = "vector",
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

## Show
if(!isGeneric("show")){
  setGeneric(name = "show",
             def = function(object) standardGeneric("show"))
}

setMethod("show", "simulation",
          function(object){show.simulation(object)})


#' Anisotropic simulation class
#'
#' @description
#' S4 class for anisotropic 2D particle movement simulation.
#'
#' @slot sz vector. Frame size of the intensity profile, number of pixels
#' contained in each frame equals \code{sz[1]} by \code{sz[2]}.
#' @slot len_t integer. Number of time steps.
#' @slot noise character. Background noise, options from ('uniform','gaussian').
#' @slot model_name character. Simulated stochastic process, options from ('BM','OU','FBM','OU+FBM').
#' @slot M integer. Number of particles.
#' @slot pxsz numeric.  Size of one pixel in unit of micron, 1 for simulated data.
#' @slot mindt numeric. Minimum lag time, 1 for simulated data.
#' @slot pos matrix. Position matrix for particle trajectory, see 'Details'.
#' @slot intensity matrix. Filled intensity profile, see 'Details'.
#' @slot num_msd matrix. Numerical mean squared displacement (MSD).
#' @slot param matrix. Parameters used to construct MSD.
#' @slot theor_msd matrix. Theoretical MSD.
#' @slot sigma_2_0 vector. Variance of background noise.
#'
#' @method show aniso_simulation
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
methods::setClass("aniso_simulation", representation(
  sz = "vector",
  len_t = "integer",
  noise = "character",
  model_name = "character",
  M = "integer",
  pxsz = "numeric",
  mindt = "numeric",
  pos = "matrix", #first M equals pos0
  intensity = "matrix",
  num_msd = "matrix",
  param = "matrix",
  theor_msd = "matrix",
  sigma_2_0 = "vector"
)
)

## Show
if(!isGeneric("show")){
  setGeneric(name = "show",
             def = function(object) standardGeneric("show"))
}

setMethod("show", "aniso_simulation",
          function(object){show.aniso_simulation(object)})
