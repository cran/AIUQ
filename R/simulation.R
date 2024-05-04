#' Simulate 2D particle movement
#'
#' @description
#' Simulate 2D particle movement from a user selected stochastic process, and
#' output intensity profiles.
#'
#' @param sz frame size of simulated image with default \code{c(200,200)}.
#' @param len_t number of time steps with default 200.
#' @param M number of particles with default 50.
#' @param model_name stochastic process simulated, options from
#' ('BM','OU','FBM','OU+FBM'), with default 'BM'.
#' @param noise background noise, options from ('uniform','gaussian'),
#' with default 'gaussian'.
#' @param I0 background intensity, value between 0 and 255, with default 20.
#' @param Imax maximum intensity at the center of the particle, value between 0
#' and 255, with default 255.
#' @param pos0 initial position for M particles, matrix with dimension M by 2.
#' @param rho correlation between successive step and previous step in O-U
#' process, value between 0 and 1, with default 0.95.
#' @param H Hurst parameter of fractional Brownian Motion, value between 0 and 1,
#' with default 0.3.
#' @param sigma_p radius of the spherical particle (3sigma_p), with default 2.
#' @param sigma_bm distance moved per time step in Brownian Motion, with default 1.
#' @param sigma_ou distance moved per time step in Ornstein–Uhlenbeck process,
#' with default 2.
#' @param sigma_fbm distance moved per time step in fractional Brownian Motion,
#' with default 2.
#' @return Returns an S4 object of class \code{simulation}.
#' @export
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
#' @examples
#' library(AIUQ)
#' # -------------------------------------------------
#' # Example 1: Simple diffusion for 200 images with
#' #            200 by 200 pixels and 50 particles
#' # -------------------------------------------------
#' sim_bm = simulation()
#' show(sim_bm)
#'
#' # -------------------------------------------------
#' # Example 2: Simple diffusion for 100 images with
#' #            100 by 100 pixels and slower speed
#' # -------------------------------------------------
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#'
#' # -------------------------------------------------
#' # Example 3: Ornstein-Uhlenbeck process
#' # -------------------------------------------------
#' sim_ou = simulation(model_name="OU")
#' show(sim_ou)
simulation <- function(sz=c(200,200), len_t=200, M=50, model_name="BM",noise="gaussian",
                       I0=20, Imax=255, pos0=matrix(NaN,nrow=M,ncol=2),rho=0.95,
                       H=0.3, sigma_p=2, sigma_bm=1,sigma_ou=2, sigma_fbm=2){

  model <- methods::new("simulation")

  #check
  len_t = as.integer(len_t)
  M = as.integer(M)
  if(length(sz)==1){
    sz=c(sz,sz)
  }
  if(length(sz)>2){
    stop("Frame size of simulated image should be a vector with length 2. \n")
  }
  if(!is.character(model_name)){
    stop("Type of stochastic process should be a character value. \n")
  }
  if(model_name!="BM" && model_name!="OU" && model_name!="FBM" && model_name!="OU+FBM"){
    stop("Type of stochastic process should be one of the type listed in help page. \n")
  }
  if(!is.character(noise)){
    stop("Type of background noise should be a character value. \n")
  }
  if(noise!="gaussian" && noise!="uniform"){
    stop("Type of background noise should be one of the type listed in help page. \n")
  }
  if(!is.numeric(I0)){
    stop("Background intensity should have numeric value. \n")
  }
  if(I0<0 || I0>255){
    stop("Background intensity should have value between 0 and 255. \n")
  }
  if(!is.numeric(Imax)){
    stop("Maximum intensity at the center of the particle should be a numeric value. \n")
  }
  if(Imax<0 || Imax>255){
    stop("Maximum intensity at the center of the particle should have value between 0 and 255. \n")
  }
  if(!is.numeric(pos0)){
    stop("Initial position for particles should be all numeric. \n")
  }
  if(nrow(pos0)!=M || ncol(pos0)!=2){
    stop("Dimension of particle initial position matrix should match M by 2. \n")
  }
  if(!is.numeric(rho)){
    stop("Correlation between steps in O-U process should be numeric. \n")
  }
  if(!is.numeric(H)){
    stop("Hurst parameter of fractional Brownian Motion should be numeric. \n")
  }
  if(H<0 || H>1){
    stop("Hurst parameter of fractional Brownian Motion should have value between 0 and 1. \n")
  }
  if(!is.numeric(sigma_p)){
    stop("Radius of the spherical particle should be numeric. \n")
  }
  if(!is.numeric(sigma_bm)){
    stop("Distance moved per time step in Brownian Motion should be numeric. \n")
  }
  if(!is.numeric(sigma_ou)){
    stop("Distance moved per time step in Ornstein Uhlenbeck process should be numeric. \n")
  }
  if(!is.numeric(sigma_fbm)){
    stop("Distance moved per time step in fractional Brownian Motion should be numeric. \n")
  }

  # Simulation particle trajectory for isotropic process
  if(sum(is.na(pos0))>=1){
    pos0 = matrix(c(sz[2]/8+0.75*sz[2]*stats::runif(M),
                    sz[1]/8+0.75*sz[1]*stats::runif(M)),nrow=M,ncol=2)
  }
  if(model_name == "BM"){
    pos = bm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma_bm)
    model@param = c(sigma_bm)
  }else if(model_name == "OU"){
    pos = ou_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma_ou,rho=rho)
    model@param = c(rho,sigma_ou)
  }else if(model_name == "FBM"){
    pos = fbm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma_fbm,H=H)
    model@param = c(sigma_fbm,H)
  }else if(model_name == "OU+FBM"){
    pos = fbm_ou_particle_intensity(pos0=pos0,M=M,len_t=len_t,H=H, rho=rho,
                                    sigma_ou = sigma_ou, sigma_fbm = sigma_fbm)
    model@param = c(rho,sigma_ou,sigma_fbm,H)
  }

  model_param = get_true_param_sim(param_truth=model@param,model_name=model_name)
  model@theor_msd = get_MSD(theta = model_param ,d_input=0:(len_t-1),model_name=model_name)
  # Fill intensity
  if(length(I0) == len_t){
    if(noise == "uniform"){
      I = matrix(stats::runif(sz[1]*sz[2]*len_t)-0.5, nrow=len_t,ncol = sz[1]*sz[2])
      I = I*I0
      model@sigma_2_0 = I0^2/12
    }else if(noise == "gaussian"){
      I = matrix(stats::rnorm(sz[1]*sz[2]*len_t), nrow=len_t,ncol = sz[1]*sz[2])
      I = I*sqrt(I0)
      model@sigma_2_0 = I0
    }
  }else if(length(I0) == 1){
    if(noise == "uniform"){
      I = matrix(I0*(stats::runif(sz[1]*sz[2]*len_t)-0.5), nrow=len_t,ncol = sz[1]*sz[2])
      model@sigma_2_0 = c(I0^2/12)
    } else if(noise == "gaussian"){
      I = matrix(sqrt(I0)*stats::rnorm(sz[1]*sz[2]*len_t), nrow=len_t,ncol = sz[1]*sz[2])
      model@sigma_2_0 = c(I0)
    }
  }

  if(length(Imax)==1){
    Ic = rep(Imax,M)
    model@intensity = fill_intensity(len_t=len_t,M=M,I=I,pos=pos,Ic=Ic,sz=sz, sigma_p=sigma_p)
  }

  model@sz = sz
  model@pxsz = 1
  model@mindt = 1
  model@len_t = len_t
  model@noise = noise
  model@M = M
  model@model_name = model_name
  model@pos = pos
  model@num_msd = numerical_msd(pos=model@pos,M=model@M,len_t=model@len_t)

  return(model)
}
