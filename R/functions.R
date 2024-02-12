#' Transform Cartesian coordinates to polar coordinates
#' @description
#' Transform ordered pair (x,y), where x and y denotes the
#' directed distance between the point and each of two
#' perpendicular lines, the x-axis and the y-axis, to polar
#' coordinate. Input x and y must have the same length.
#'
#'
#' @param x a vector of x-coordinates
#' @param y a vector of y-coordinates
#'
#' @return A data frame with 2 variables, where r is the
#' directed distance from a point designed as the pole, and
#' theta represents the angle, in radians, between the pole and the point.
#' @export
#' @author \packageAuthor{AIUQ}
#' @concept cart2pol
#' @examples
#' library(AIUQ)
#'
#' # Input in Cartesian coordinates
#' (x <- rep(1:3,each = 3))
#' (y <- rep(1:3,3))
#'
#' # Data frame with polar coordinates
#' (polar <- cart2polar(x, y))
#'
#' @keywords internal
cart2polar <- function(x, y) {
  if(length(x)!=length(y)){
    stop("Length of points in each coordinates shoule be the same. \n")
  }
  if(!is.numeric(x) || !is.numeric(y)){
    stop("Inputs shoule have numeric value. \n")
  }
  data.frame(theta = atan2(y, x),r = sqrt(x^2 + y^2))
}

#' Transform intensity profile into SS_T matrix
#'
#' @description
#' Transform intensity profile with different formats, ('SST_array','T_SS_mat',
#' 'SS_T_mat','S_ST_mat'), space by space by time array, time by (space by space) matrix,
#' (space by space) by time matrix, or space by (space by time) matrix, into
#' 'SS_T_mat'. In addition, crop each frame with odd frame size.
#'
#' @param intensity intensity profile, array or matrix
#' @param intensity_str structure of the original intensity
#' profile, options from ('SST_array','T_SS_mat','SS_T_mat','S_ST_mat')
#' @param square a logical evaluating to TRUE or FALSE indicating whether or not
#' to crop each frame into a square such that frame size in x direction equals
#' frame size in y direction with \code{sz_x=sz_y}
#' @param sz frame size of each frame. A vector of length 2 with frame size in
#' y/(row) direction, and frame size in x/(column) direction, with default \code{NA}.
#'
#' @return A matrix of transformed intensity profile.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' # -------------------------------------------------
#' # Example 1: Transform T_SS_mat into SS_T_mat, each
#' #             frame contains number 1-9
#' # -------------------------------------------------
#' (m <- matrix(rep(1:9,4),4,9,byrow=TRUE))
#' intensity_format_transform(m,intensity_str="T_SS_mat",sz=c(4,9))
#'
#' # -------------------------------------------------
#' # Example 2: Transform SST_array into SS_T_mat, each
#' #             frame contains number 1-9
#' # -------------------------------------------------
#' (m <- array(rep(1:9,4),dim=c(3,3,4)))
#' intensity_format_transform(m,intensity_str="SST_array")
#'
#' @keywords internal
intensity_format_transform<-function(intensity,intensity_str, square=FALSE, sz=NA){
  if(intensity_str=='SST_array'){#most real experimental structure
    if(square==TRUE){
      sz_x = sz_y = min(dim(intensity)[1],dim(intensity)[2])
      if(sz_x%%2 == 0){ #if even frame size, crop into odd
        sz_x=sz_x-1
        sz_y=sz_y-1
      }
    }else{
      sz_y = dim(intensity)[1]
      sz_x = dim(intensity)[2]
      if(sz_y%%2==0){
        sz_y = sz_y-1
      }
      if(sz_x%%2==0){
        sz_x = sz_x-1
      }
    }
    len_t = dim(intensity)[3]
    intensity_transform = matrix(NA,sz_x*sz_y,len_t)
    for(i in 1:len_t){
      intensity_transform[,i] = as.numeric(intensity[1:sz_y,1:sz_x,i])
    }

  }else if(intensity_str=='T_SS_mat'){#Simulated data using AIUQ simulation class
    if(sum(is.na(sz))>0){
      stop("Please give a vector of sz that contains frame size of each frame.")
    }
    intensity=as.matrix(intensity)
    if(square==TRUE){
      sz_x=sz_y=min(sz)
      if(sz_x%%2==0){
        sz_x=sz_x-1
        sz_y=sz_y-1
      }
    }else{
      sz_y = sz[1]
      sz_x = sz[2]
      if(sz_y%%2==0){
        sz_y = sz_y-1
      }
      if(sz_x%%2==0){
        sz_x = sz_x-1
      }
    }
    len_t = dim(intensity)[1]
    intensity_transform=matrix(NA,sz_x*sz_y,len_t)
    for(i in 1:len_t){
      intensity_mat=matrix(intensity[i,],sz[1],sz[2])
      intensity_transform[,i]=as.vector(intensity_mat[1:sz_y,1:sz_x])
    }

  }else if(intensity_str=='S_ST_mat'){
    if(sum(is.na(sz))>0){
      stop("Please give a vector of sz that contains frame size of each frame.")
    }
    intensity=as.matrix(intensity)
    if(square==TRUE){
      sz_x=sz_y=min(sz)
      if(sz_x%%2==0){
        sz_x=sz_x-1
        sz_y=sz_y-1
      }
    }else{
      sz_y = sz[1]
      sz_x = sz[2]
      if(sz_y%%2==0){
        sz_y=sz_y-1
      }
      if(sz_x%%2==0){
        sz_x=sz_x-1
      }
    }
    len_t = dim(intensity)[2]/sz[2]
    intensity_transform=matrix(NA,sz_x*sz_y,len_t)
    for(i in 1:len_t){
      selected_x = ((1+sz[2]*(i-1)):(sz[2]+sz[2]*(i-1)))[1:sz_x]
      intensity_transform[,i]=
        as.vector(intensity[1:sz_y,selected_x])
    }

  }else if(intensity_str=='SS_T_mat'){
    if(sum(is.na(sz))>0){
      stop("Please give a vector of sz that contains frame size of each frame.")
    }
    intensity=as.matrix(intensity)
    if(square==TRUE){
      sz_x=sz_y=min(sz)
      if(sz_x%%2==0){
        sz_x=sz_x-1
        sz_y=sz_y-1
      }
    }else{
      sz_y = sz[1]
      sz_x = sz[2]
      if(sz_y%%2==0){
        sz_y = sz_y-1
      }
      if(sz_x%%2==0){
        sz_x=sz_x-1
      }
    }
    len_t = dim(intensity)[2]
    intensity_transform=matrix(NA,sz_x*sz_y,len_t)
    for(i in 1:len_t){
      intensity_mat=matrix(intensity[,i],sz[1],sz[2])
      intensity_transform[,i]=as.vector(intensity_mat[1:sz_y,1:sz_x])
    }

  }
  intensity_transform_list = list()
  intensity_transform_list$sz_x = sz_x
  intensity_transform_list$sz_y = sz_y
  intensity_transform_list$intensity = intensity_transform
  return(intensity_transform_list)
}

#' 2D Fourier transformation and calculate wave number
#' @description
#' Perform 2D fast Fourier transformation on SS_T matrix, record frame size for
#' each frame, total number of frames, and sequence of lag times. Calculate and
#' record circular wave number for complete q ring.
#'
#' @param intensity_list intensity profile, SS_T matrix
#' @param pxsz size of one pixel in unit of micron
#' @param mindt minimum lag time
#'
#' @return A list object containing transformed intensity profile in reciprocal
#' space and corresponding parameters.
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
#' @keywords internal
FFT2D<-function(intensity_list,pxsz,mindt){
  sz_x = intensity_list$sz_x
  sz_y = intensity_list$sz_y
  intensity = intensity_list$intensity
  len_t = dim(intensity)[2]
  I_q_matrix = matrix(NA,sz_x*sz_y,len_t)
  for(i in 1:len_t){
    I_q_matrix[,i] = as.vector(fftwtools::fftw2d(matrix(intensity[,i],sz_y,sz_x)))
  }


  ans_list = list()
  #ans_list$sz = sz
  ans_list$sz_x = sz_x
  ans_list$sz_y = sz_y
  ans_list$len_q = length(1:((max(sz_x,sz_y)-1)/2))
  ans_list$len_t = len_t
  ans_list$I_q_matrix = I_q_matrix
  ans_list$q = (1:((max(sz_x,sz_y)-1)/2))*2*pi/(max(sz_x,sz_y)*pxsz)
  ans_list$input = mindt*(1:(len_t))
  ans_list$d_input = ans_list$input[1:length(ans_list$input)]-ans_list$input[1] ##delta t, including zero

  return(ans_list)
}

#' fftshift
#'
#' @description
#' Rearranges a 2D Fourier transform x by shifting the zero-frequency component
#' to the center of the matrix.
#'
#'
#' @param x square matrix input with odd number of rows and columns
#' @param dim shift method. See 'Details'.
#'
#' @return Shifted matrix.
#'
#' @details
#' By default, \code{dim=-1}, swaps the first quadrant of x with the
#' third, and the second quadrant with the fourth. If \code{dim=1}, swaps rows 1
#' to middle with rows (middle+1) to end. If \code{dim=2}, swaps columns 1
#' to middle with columns (middle+1) to end. If \code{dim=3}, reverse fftshift.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#'
#' (m <- matrix(0:8,3,3))
#' fftshift(m)
#'
#' @keywords internal
fftshift <- function(x, dim = -1) {

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  swap_up_down <- function(x) {
    rows_half <- ceiling(rows/2)
    return(rbind(x[((rows_half+1):rows), (1:cols)], x[(1:rows_half), (1:cols)]))
  }

  swap_left_right <- function(x) {
    cols_half <- ceiling(cols/2)
    return(cbind(x[1:rows, ((cols_half+1):cols)], x[1:rows, 1:cols_half]))
  }

  swap_up_down_reverse <- function(x) {
    rows_half <- ceiling(rows/2)
    return(rbind(x[((rows_half):rows), (1:cols)], x[1:(rows_half-1), (1:cols)]))
  }

  swap_left_right_reverse <- function(x) {
    cols_half <- ceiling(cols/2)
    return(cbind(x[1:rows, ((cols_half):cols)], x[1:rows, 1:(cols_half-1)]))
  }


  if (dim == -1) {
    x <- swap_up_down(x)
    return(swap_left_right(x))
  }
  else if (dim == 1) {
    return(swap_up_down(x))
  }
  else if (dim == 2) {
    return(swap_left_right(x))
  }else if(dim == 3){
    x <- swap_up_down_reverse(x)
    return(swap_left_right_reverse(x))

  }
  else {
    stop("Invalid dimension parameter")
  }
}

#' Construct MSD and MSD gradient with transformed parameters
#'
#' @description
#' Construct mean squared displacement (MSD) and its gradient for a given
#' stochastic process or a user defined MSD and gradient structure.
#'
#' @param theta transformed parameters in MSD function for MLE estimation
#' @param d_input sequence of lag times
#' @param model_name model name for the process, options from  ('BM','OU','FBM',
#' 'OU+FBM', 'user_defined').
#' @param msd_fn user defined MSD structure, a function of \code{theta} and \code{d_input}
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{theta} and \code{d_input}
#'
#' @return A list of two variables, MSD and MSD gradient.
#' @details
#' Note for non \code{user_defined} model, \code{msd_fn} and \code{msd_grad_fn}
#' are not needed. For Brownian Motion, the MSD follows
#' \deqn{MSD_{BM}(\Delta t) = \theta_1\Delta t= 4D\Delta t}{%MSD_{BM}(\Delta t) = \theta_1\Delta t= 4D\Delta t}
#' where \code{D} is the diffusion coefficient.
#'
#' For Ornstein–Uhlenbeck process,  the MSD follows
#' \deqn{MSD_{OU}(\Delta t) = \theta_2(1-\frac{\theta_1}{1+\theta_1}^{\Delta t})}{%MSD_{OU}(\Delta t) = \theta_2(1-\frac{\theta_1}{1+\theta_1}^{\Delta t})}
#' where \eqn{\frac{\theta_1}{1+\theta_1}=\rho}{%\frac{\theta_1}{1+\theta_1}=\rho}
#' is the correlation with previous steps.
#'
#' For fractional Brownian Motion,  the MSD follows
#' \deqn{MSD_{FBM}(\Delta t) =\theta_1\Delta t^{\frac{2\theta_2}{1+\theta_2}}}{%MSD_{FBM}(\Delta t) =\theta_1\Delta t^{\frac{2\theta_2}{1+\theta_2}}}
#' where \eqn{\frac{2\theta_2}{1+\theta_2}=2H}{%\frac{2\theta_2}{1+\theta_2}=2H}
#' with \code{H} is the the Hurst parameter.
#'
#' For 'OU+FBM', the MSD follows
#' \deqn{MSD_{OU+FBM}(\Delta t) = \theta_2(1-\frac{\theta_1}{1+\theta_1}^{\Delta t})+\theta_3\Delta t^{\frac{2\theta_4}{1+\theta_4}}}{%MSD_{OU+FBM}(\Delta t) = \theta_2(1-\frac{\theta_1}{1+\theta_1}^{\Delta t})+\theta_3\Delta t^{\frac{2\theta_4}{1+\theta_4}}}
#'
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
#' @examples
#' library(AIUQ)
#' msd_fn <- function(param, d_input){
#'   beta = 2*param[1]^2
#'   MSD = beta*d_input
#' }
#' msd_grad_fn <- function(param, d_input){
#'   MSD_grad = 4*param[1]*d_input
#' }
#'
#' theta = 2
#' d_input = 0:10
#' model_name = "user_defined"
#'
#' MSD_list = get_MSD_with_grad(theta=theta,d_input=d_input,
#'                              model_name=model_name,msd_fn=msd_fn,
#'                              msd_grad_fn=msd_grad_fn)
#' MSD_list$msd
#' MSD_list$msd_grad
#'
#' @keywords internal
get_MSD_with_grad<-function(theta,d_input,model_name,msd_fn=NA,msd_grad_fn=NA){
  if(model_name=='user_defined' || model_name=='user_defined_anisotropic'){
    if(is.function(msd_grad_fn)==T){
      MSD = msd_fn(theta, d_input)
      MSD_grad = msd_grad_fn(theta, d_input)
    }else{
      MSD = msd_fn(theta, d_input)
      MSD_grad = NA
    }
  }else if(model_name=='BM'||model_name=='BM_anisotropic'){
    beta = theta[1]
    MSD = beta*d_input
    MSD_grad = as.matrix(d_input)
  }else if(model_name=='FBM'||model_name=='FBM_anisotropic'){
    beta = theta[1]
    alpha = 2*theta[2]/(1+theta[2])
    MSD = beta*d_input^alpha
    MSD_grad = cbind(d_input^alpha, beta*c(0,log(d_input[-1]))*(d_input^alpha))
  }else if(model_name=='OU'||model_name=='OU_anisotropic'){
    rho = theta[1]/(1+theta[1])
    amplitude = theta[2]
    MSD = (amplitude*(1-rho^d_input))
    MSD_grad = cbind(-amplitude*d_input*(rho^(d_input-1)), (1-rho^d_input))
  }else if(model_name=='OU+FBM'||model_name=='OU+FBM_anisotropic'){
    rho = theta[1]/(1+theta[1])
    amplitude = theta[2]
    beta = theta[3]
    alpha = 2*theta[4]/(1+theta[4])
    MSD = beta*d_input^alpha+(amplitude*(1-rho^d_input))
    MSD_grad = cbind(-amplitude*d_input*(rho^(d_input-1)),
                     (1-rho^d_input),
                     d_input^alpha,beta*c(0,log(d_input[-1]))*(d_input^alpha))
  }else if(model_name=='VFBM'||model_name=='VFBM_anisotropic'){
    a = theta[1]
    b = theta[2]
    c = theta[3]/(1+theta[3])
    d = theta[4]/(1+theta[4])
    MSD = a*d_input^((c*d_input)/(1+b*d_input)+d)
    MSD_grad = cbind(d_input^((c*d_input)/(1+b*d_input)+d),
                     -a*d_input^((c*d_input)/(1+b*d_input)+d)*c(0,log(d_input[-1]))*c*(1+b*d_input)^(-2)*d_input^2,
                     a*d_input^((c*d_input)/(1+b*d_input)+d)*c(0,log(d_input[-1]))*d_input/(1+b*d_input),
                     a*d_input^((c*d_input)/(1+b*d_input)+d)*c(0,log(d_input[-1])))
  }
  msd_list = list()
  msd_list$msd = MSD
  msd_list$msd_grad = MSD_grad
  return(msd_list)
}

#' Transform parameters in simulation class to parameters structure in MSD function
#' @description
#' Transform parameters in \code{simulation} class to parameters contained in MSD
#' function with structure \code{theta} in \code{\link{get_MSD}}. Prepare for
#' truth MSD construction.
#'
#' @param param_truth parameters used in \code{simulation} class
#' @param model_name stochastic process used in \code{simulation}, options from
#' ('BM','OU','FBM','OU+FBM')
#'
#' @return A vector of parameters contained in MSD with structure \code{theta} in
#' \code{\link{get_MSD}}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' # Simulate simple diffusion for 100 images with 100 by 100 pixels and
#' # distance moved per time step is 0.5
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#' get_true_param_sim(param_truth=sim_bm@param,model_name=sim_bm@model_name)
#'
#' @keywords internal
get_true_param_sim<-function(param_truth,model_name){
  if(model_name=='BM'){
    beta = 2*param_truth[1]^2
    param = c(beta)
  }else if(model_name=='FBM'){
    beta = 2*param_truth[1]^2
    alpha = 2*param_truth[2]
    param = c(beta, alpha)
  }else if(model_name=='OU'){
    rho = param_truth[1]
    amplitude = 4*param_truth[2]^2
    param = c(rho,amplitude)
  }else if(model_name=='OU+FBM'){
    rho = param_truth[1]
    amplitude = 4*param_truth[2]^2
    beta = 2*param_truth[1]^2
    alpha = 2*param_truth[2]
    param = c(rho,amplitude,beta,alpha)
  }
  return(param)
}

#' Transform parameters in anisotropic simulation class to parameters structure in MSD function
#' @description
#' Transform parameters in \code{aniso_simulation} class to parameters contained in MSD
#' function with structure \code{theta} in \code{\link{get_MSD}}. Prepare for
#' truth MSD construction.
#'
#' @param param_truth parameters used in \code{aniso_simulation} class
#' @param model_name stochastic process used in \code{aniso_simulation}, options from
#' ('BM','OU','FBM','OU+FBM')
#'
#' @return A vector of parameters contained in MSD with structure \code{theta} in
#' \code{\link{get_MSD}}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' # Simulate simple diffusion for 100 images with 100 by 100 pixels and
#' # distance moved per time step is 0.5
#' aniso_sim_bm = aniso_simulation(sz=100,len_t=100,sigma_bm=c(0.5,0.1))
#' show(aniso_sim_bm)
#' get_true_param_aniso_sim(param_truth=aniso_sim_bm@param,model_name=aniso_sim_bm@model_name)
#'
#' @keywords internal
get_true_param_aniso_sim<-function(param_truth,model_name){
  if(model_name=='BM'){
    beta = param_truth[1]^2
    param = c(beta)
  }else if(model_name=='FBM'){
    beta = param_truth[1]^2
    alpha = 2*param_truth[2]
    param = c(beta, alpha)
  }else if(model_name=='OU'){
    rho = param_truth[1]
    amplitude = 2*param_truth[2]^2
    param = c(rho,amplitude)
  }else if(model_name=='OU+FBM'){
    rho = param_truth[1]
    amplitude = 2*param_truth[2]^2
    beta = param_truth[3]^2
    alpha = 2*param_truth[4]
    param = c(rho,amplitude,beta,alpha)
  }
  return(param)
}


#' Transform parameters estimated in SAM class to parameters structure in MSD function
#' @description
#' Transform parameters estimated using Maximum Likelihood Estimation (MLE) in
#' \code{SAM} class, to parameters contained in MSD with structure \code{theta}
#' in \code{\link{get_MSD}}.
#'
#' @param theta estimated parameters through MLE
#' @param model_name fitted stochastic process, options from ('BM','OU','FBM','OU+FBM')
#'
#' @return A vector of estimated parameters after transformation with structure
#' \code{theta} in \code{\link{get_MSD}}.
#'
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
#' @keywords internal
get_est_param<-function(theta,model_name){
  if(model_name=='BM'){
    beta = theta[1]
    est_param = c(beta)
  }else if(model_name=='FBM'){
    beta = theta[1]
    alpha = 2*theta[2]/(1+theta[2])
    est_param = c(beta, alpha)
  }else if(model_name=='OU'){
    rho = theta[1]/(1+theta[1])
    amplitude = theta[2]
    est_param = c(rho,amplitude)
  }else if(model_name=='OU+FBM'){
    rho = theta[1]/(1+theta[1])
    amplitude = theta[2]
    beta = theta[3]
    alpha = 2*theta[4]/(1+theta[4])
    est_param = c(rho,amplitude,beta,alpha)
  }else if(model_name=='user_defined'){
    est_param = theta
  }else if(model_name=='VFBM'){
    a = theta[1]
    b = theta[2]
    c = theta[3]/(1+theta[3])
    d = theta[4]/(1+theta[4])
    est_param = c(a,b,c,d)
  }
  return(est_param)
}


#' Construct MSD
#' @description
#' Construct estimated mean squared displacement (MSD) for a given stochastic process.
#'
#' @param theta parameters in MSD function
#' @param d_input sequence of lag times
#' @param model_name model name for the process, options from ('BM','OU','FBM',
#' 'OU+FBM','user_defined')
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#'
#' @return A vector of MSD values for a given sequence of lag times.
#' @details
#' For Brownian Motion, the MSD follows
#' \deqn{MSD_{BM}(\Delta t) = \theta_1\Delta t= 4D\Delta t}{%MSD_{BM}(\Delta t) = \theta_1\Delta t= 4D\Delta t}
#' where \code{D} is the diffusion coefficient.
#'
#' For Ornstein–Uhlenbeck process,  the MSD follows
#' \deqn{MSD_{OU}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})}{%MSD_{OU}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})}
#' where \eqn{\theta_1=\rho}{%\theta_1=\rho}
#' is the correlation with previous steps.
#'
#' For fractional Brownian Motion,  the MSD follows
#' \deqn{MSD_{FBM}(\Delta t) =\theta_1\Delta t^{\theta_2}}{%MSD_{FBM}(\Delta t) =\theta_1\Delta t^{\theta_2}}
#' where \eqn{\theta_2=2H}{%\theta_2=2H} with \code{H} is the the Hurst parameter.
#'
#' For 'OU+FBM', the MSD follows
#' \deqn{MSD_{OU+FBM}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})+\theta_3\Delta t^{\theta_4}}{%MSD_{OU+FBM}(\Delta t) = \theta_2(1-\theta_1^{\Delta t})+\theta_3\Delta t^{\theta_4}}
#'
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
#' @examples
#' library(AIUQ)
#' # Construct MSD for BM
#' get_MSD(theta=0.2,d_input=0:100,model_name='BM')
#'
#' @keywords internal
get_MSD<-function(theta,d_input,model_name, msd_fn=NA){
  if(model_name=='BM'){
    beta = theta[1]
    MSD = beta*d_input
  }else if(model_name=='FBM'){
    beta = theta[1]
    alpha = theta[2]
    MSD = beta*d_input^alpha
  }else if(model_name=='OU'){
    rho = theta[1]
    amplitude = theta[2]
    MSD = amplitude*(1-rho^d_input)
  }else if(model_name=='OU+FBM'){
    rho = theta[1]
    amplitude = theta[2]
    beta = theta[3]
    alpha = theta[4]
    MSD = beta*d_input^alpha+(amplitude*(1-rho^d_input))
  }else if(model_name=='user_defined'){
    MSD = msd_fn(theta, d_input)
  }else if(model_name=='VFBM'){
    a = theta[1]
    b = theta[2]
    c = theta[3]
    d = theta[4]
    power = (c*d_input)/(b*d_input+1)+d
    MSD = a*d_input^power
  }
  return(MSD)
}

#' Construct 95% confidence interval
#' @description
#' This function construct the lower and upper bound for 95% confidence interval
#' of estimated parameters and mean squared displacement(MSD) for a given model.
#' See 'References'.
#'
#' @param param_uq_range lower and upper bound for natural logorithm of
#' parameters in the fitted model using \code{AIUQ} method in \code{SAM} class
#' @param model_name model for constructing MSD, options from ('BM','OU',
#' 'FBM','OU+FBM', 'user_defined')
#' @param d_input sequence of lag times
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#'
#' @return A list of lower and upper bound for 95% confidence interval
#' of estimated parameters and MSD for a given model.
#'
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
#' @keywords internal
get_est_parameters_MSD_SAM_interval <- function(param_uq_range,model_name,d_input,msd_fn=NA){

  theta=exp(param_uq_range[,-dim(param_uq_range)[2]])
  sigma_2_0_est=exp(param_uq_range[,dim(param_uq_range)[2]])
  sigma_2_0_est_lower=sigma_2_0_est[1]
  sigma_2_0_est_upper=sigma_2_0_est[2]

  est_parameters=NA;

  if(model_name=='BM'){
    beta_lower=theta[1] ##only 1 param
    beta_upper=theta[2]

    MSD_lower=beta_lower*d_input
    MSD_upper=beta_upper*d_input

    est_parameters_lower=c(beta_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(beta_upper,sigma_2_0_est_upper);

  }else if(model_name=='FBM'){
    beta_lower=theta[1,1]
    beta_upper=theta[2,1]

    alpha_lower=2*theta[1,2]/(1+theta[1,2])
    alpha_upper=2*theta[2,2]/(1+theta[2,2])

    MSD_lower=rep(NA,length(d_input))
    MSD_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)
    if(length(index_less_than_1)>0){
      MSD_lower[index_less_than_1]=beta_lower*d_input[index_less_than_1]^{alpha_upper}
      MSD_upper[index_less_than_1]=beta_upper*d_input[index_less_than_1]^{alpha_lower}
      MSD_lower[-index_less_than_1]=beta_lower*d_input[-index_less_than_1]^{alpha_lower}
      MSD_upper[-index_less_than_1]=beta_upper*d_input[-index_less_than_1]^{alpha_upper}

    }else{
      MSD_lower=beta_lower*d_input^{alpha_lower}
      MSD_upper=beta_upper*d_input^{alpha_upper}

    }
    est_parameters_lower=c(beta_lower,alpha_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(beta_upper,alpha_upper,sigma_2_0_est_upper);
  }else if(model_name=='OU'){ #seems okay
    rho_lower=theta[1,1]/(1+theta[1,1])
    rho_upper=theta[2,1]/(1+theta[2,1])

    amplitude_lower=theta[1,2]
    amplitude_upper=theta[2,2]


    MSD_lower=(amplitude_lower*(1-rho_lower^d_input))
    MSD_upper=(amplitude_upper*(1-rho_upper^d_input))

    est_parameters_lower=c(rho_lower,amplitude_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(rho_upper,amplitude_upper,sigma_2_0_est_upper);

  }else if(model_name=='OU+FBM'){
    rho_lower=theta[1,1]/(1+theta[1,1])
    rho_upper=theta[2,1]/(1+theta[2,1])

    amplitude_lower=theta[1,2]
    amplitude_upper=theta[2,2]

    beta_lower=theta[1,3]
    beta_upper=theta[2,3]

    alpha_lower=2*theta[1,4]/(1+theta[1,4])
    alpha_upper=2*theta[2,4]/(1+theta[2,4])

    ####change of uq, need to test
    MSD_lower=rep(NA,length(d_input))
    MSD_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)
    if(length(index_less_than_1)>0){
      MSD_lower[index_less_than_1]=beta_lower*d_input[index_less_than_1]^{alpha_upper}
      MSD_upper[index_less_than_1]=beta_upper*d_input[index_less_than_1]^{alpha_lower}
      MSD_lower[-index_less_than_1]=beta_lower*d_input[-index_less_than_1]^{alpha_lower}
      MSD_upper[-index_less_than_1]=beta_upper*d_input[-index_less_than_1]^{alpha_upper}

    }else{
      MSD_lower=beta_lower*d_input^{alpha_lower}
      MSD_upper=beta_upper*d_input^{alpha_upper}

    }

    MSD_lower=MSD_lower+(amplitude_lower*(1-rho_lower^d_input))
    MSD_upper=MSD_upper+(amplitude_upper*(1-rho_upper^d_input))

    est_parameters_lower=c(rho_lower,amplitude_lower,beta_lower,alpha_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(rho_upper,amplitude_upper,beta_upper,alpha_upper,sigma_2_0_est_upper);

  }else if(model_name=='user_defined'){
    if(is.matrix(theta)){
      theta_lower=theta[1,]
      theta_upper=theta[2,]


      MSD_lower=msd_fn(theta_lower,d_input)
      MSD_upper=msd_fn(theta_upper,d_input)

      est_parameters_lower=c(theta_lower,sigma_2_0_est_lower)
      est_parameters_upper=c(theta_upper,sigma_2_0_est_upper)
    }else{
      theta_lower=theta[1]
      theta_upper=theta[2]


      MSD_lower=msd_fn(theta_lower,d_input)
      MSD_upper=msd_fn(theta_upper,d_input)

      est_parameters_lower=c(theta_lower,sigma_2_0_est_lower)
      est_parameters_upper=c(theta_upper,sigma_2_0_est_upper)
    }

  }
  ans_list=list()
  ans_list$est_parameters_lower=est_parameters_lower
  ans_list$est_parameters_upper=est_parameters_upper
  ans_list$MSD_lower=MSD_lower
  ans_list$MSD_upper=MSD_upper

  return(ans_list)
}

#' Construct 95% confidence interval for anisotropic processes
#' @description
#' This function construct the lower and upper bound for 95% confidence interval
#' of estimated parameters and mean squared displacement(MSD) for a given
#' anisotropic model. See 'References'.
#'
#' @param param_uq_range lower and upper bound for natural logorithm of
#' parameters in the fitted model using \code{AIUQ} method in \code{aniso_SAM} class
#' @param model_name model for constructing MSD, options from ('BM','OU',
#' 'FBM','OU+FBM', 'user_defined')
#' @param d_input sequence of lag times
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#'
#' @return A list of lower and upper bound for 95% confidence interval
#' of estimated parameters and MSD for a given model.
#'
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
#' @keywords internal
get_est_parameters_MSD_SAM_interval_anisotropic <- function(param_uq_range,model_name,d_input,msd_fn=NA){

  theta=exp(param_uq_range[,-dim(param_uq_range)[2]])
  sigma_2_0_est=exp(param_uq_range[,dim(param_uq_range)[2]])
  sigma_2_0_est_lower=sigma_2_0_est[1]
  sigma_2_0_est_upper=sigma_2_0_est[2]

  est_parameters=NA

  if(model_name=='BM'){
    beta_x_lower=theta[1,1] ##only 1 param
    beta_x_upper=theta[2,1]
    beta_y_lower=theta[1,2] ##only 1 param
    beta_y_upper=theta[2,2]

    MSD_x_lower=beta_x_lower*d_input
    MSD_x_upper=beta_x_upper*d_input
    MSD_y_lower=beta_y_lower*d_input
    MSD_y_upper=beta_y_upper*d_input

    est_parameters_lower=c(beta_x_lower,beta_y_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(beta_x_upper,beta_y_upper,sigma_2_0_est_upper);

  }else if(model_name=='FBM'){
    beta_x_lower=theta[1,1]
    beta_x_upper=theta[2,1]

    alpha_x_lower=2*theta[1,2]/(1+theta[1,2])
    alpha_x_upper=2*theta[2,2]/(1+theta[2,2])

    MSD_x_lower=rep(NA,length(d_input))
    MSD_x_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)
    if(length(index_less_than_1)>0){
      MSD_x_lower[index_less_than_1]=beta_x_lower*d_input[index_less_than_1]^{alpha_x_upper}
      MSD_x_upper[index_less_than_1]=beta_x_upper*d_input[index_less_than_1]^{alpha_x_lower}
      MSD_x_lower[-index_less_than_1]=beta_x_lower*d_input[-index_less_than_1]^{alpha_x_lower}
      MSD_x_upper[-index_less_than_1]=beta_x_upper*d_input[-index_less_than_1]^{alpha_x_upper}

    }else{
      MSD_x_lower=beta_x_lower*d_input^{alpha_x_lower}
      MSD_x_upper=beta_x_upper*d_input^{alpha_x_upper}

    }

    beta_y_lower=theta[1,3]
    beta_y_upper=theta[2,3]

    alpha_y_lower=2*theta[1,4]/(1+theta[1,4])
    alpha_y_upper=2*theta[2,4]/(1+theta[2,4])

    MSD_y_lower=rep(NA,length(d_input))
    MSD_y_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)
    if(length(index_less_than_1)>0){
      MSD_y_lower[index_less_than_1]=beta_y_lower*d_input[index_less_than_1]^{alpha_y_upper}
      MSD_y_upper[index_less_than_1]=beta_y_upper*d_input[index_less_than_1]^{alpha_y_lower}
      MSD_y_lower[-index_less_than_1]=beta_y_lower*d_input[-index_less_than_1]^{alpha_y_lower}
      MSD_y_upper[-index_less_than_1]=beta_y_upper*d_input[-index_less_than_1]^{alpha_y_upper}

    }else{
      MSD_y_lower=beta_y_lower*d_input^{alpha_y_lower}
      MSD_y_upper=beta_y_upper*d_input^{alpha_y_upper}

    }

    est_parameters_lower=c(beta_x_lower,alpha_x_lower,beta_y_lower,alpha_y_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(beta_x_upper,alpha_x_upper,beta_y_upper,alpha_y_upper,sigma_2_0_est_upper);

  }else if(model_name=='OU'){ #seems okay
    rho_x_lower=theta[1,1]/(1+theta[1,1])
    rho_x_upper=theta[2,1]/(1+theta[2,1])

    amplitude_x_lower=theta[1,2]
    amplitude_x_upper=theta[2,2]


    MSD_x_lower=(amplitude_x_lower*(1-rho_x_lower^d_input))
    MSD_x_upper=(amplitude_x_upper*(1-rho_x_upper^d_input))

    rho_y_lower=theta[1,3]/(1+theta[1,3])
    rho_y_upper=theta[2,3]/(1+theta[2,3])

    amplitude_y_lower=theta[1,4]
    amplitude_y_upper=theta[2,4]


    MSD_y_lower=(amplitude_y_lower*(1-rho_y_lower^d_input))
    MSD_y_upper=(amplitude_y_upper*(1-rho_y_upper^d_input))

    est_parameters_lower=c(rho_x_lower,amplitude_x_lower,rho_y_lower,amplitude_y_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(rho_x_upper,amplitude_x_upper,rho_y_upper,amplitude_y_upper,sigma_2_0_est_upper);

  }else if(model_name=='OU+FBM'){
    rho_x_lower=theta[1,1]/(1+theta[1,1])
    rho_x_upper=theta[2,1]/(1+theta[2,1])

    amplitude_x_lower=theta[1,2]
    amplitude_x_upper=theta[2,2]

    beta_x_lower=theta[1,3]
    beta_x_upper=theta[2,3]

    alpha_x_lower=2*theta[1,4]/(1+theta[1,4])
    alpha_x_upper=2*theta[2,4]/(1+theta[2,4])

    ####change of uq, need to test
    MSD_x_lower=rep(NA,length(d_input))
    MSD_x_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)
    if(length(index_less_than_1)>0){
      MSD_x_lower[index_less_than_1]=beta_x_lower*d_input[index_less_than_1]^{alpha_x_upper}
      MSD_x_upper[index_less_than_1]=beta_x_upper*d_input[index_less_than_1]^{alpha_x_lower}
      MSD_x_lower[-index_less_than_1]=beta_x_lower*d_input[-index_less_than_1]^{alpha_x_lower}
      MSD_x_upper[-index_less_than_1]=beta_x_upper*d_input[-index_less_than_1]^{alpha_x_upper}

    }else{
      MSD_x_lower=beta_x_lower*d_input^{alpha_x_lower}
      MSD_x_upper=beta_x_upper*d_input^{alpha_x_upper}

    }

    MSD_x_lower=MSD_x_lower+(amplitude_x_lower*(1-rho_x_lower^d_input))
    MSD_x_upper=MSD_x_upper+(amplitude_x_upper*(1-rho_x_upper^d_input))

    rho_y_lower=theta[1,5]/(1+theta[1,5])
    rho_y_upper=theta[2,5]/(1+theta[2,5])

    amplitude_y_lower=theta[1,6]
    amplitude_y_upper=theta[2,6]

    beta_y_lower=theta[1,7]
    beta_y_upper=theta[2,7]

    alpha_y_lower=2*theta[1,8]/(1+theta[1,8])
    alpha_y_upper=2*theta[2,8]/(1+theta[2,8])

    ####change of uq, need to test
    MSD_y_lower=rep(NA,length(d_input))
    MSD_y_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)
    if(length(index_less_than_1)>0){
      MSD_y_lower[index_less_than_1]=beta_y_lower*d_input[index_less_than_1]^{alpha_y_upper}
      MSD_y_upper[index_less_than_1]=beta_y_upper*d_input[index_less_than_1]^{alpha_y_lower}
      MSD_y_lower[-index_less_than_1]=beta_y_lower*d_input[-index_less_than_1]^{alpha_y_lower}
      MSD_y_upper[-index_less_than_1]=beta_y_upper*d_input[-index_less_than_1]^{alpha_y_upper}

    }else{
      MSD_y_lower=beta_y_lower*d_input^{alpha_y_lower}
      MSD_y_upper=beta_y_upper*d_input^{alpha_y_upper}

    }

    MSD_y_lower=MSD_y_lower+(amplitude_y_lower*(1-rho_y_lower^d_input))
    MSD_y_upper=MSD_y_upper+(amplitude_y_upper*(1-rho_y_upper^d_input))

    est_parameters_lower=c(rho_x_lower,amplitude_x_lower,beta_x_lower,alpha_x_lower,
                           rho_y_lower,amplitude_y_lower,beta_y_lower,alpha_y_lower,sigma_2_0_est_lower);
    est_parameters_upper=c(rho_x_upper,amplitude_x_upper,beta_x_upper,alpha_x_upper,
                           rho_y_upper,amplitude_y_upper,beta_y_upper,alpha_y_upper,sigma_2_0_est_upper);

  }else if(model_name=='user_defined'){
      theta_x_lower=theta[1,]
      theta_x_upper=theta[2,]
      theta_y_lower=theta[3,]
      theta_y_upper=theta[4,]

      MSD_x_lower=msd_fn(theta_x_lower,d_input)
      MSD_x_upper=msd_fn(theta_x_upper,d_input)
      MSD_y_lower=msd_fn(theta_y_lower,d_input)
      MSD_y_upper=msd_fn(theta_y_upper,d_input)
      est_parameters_lower=c(theta_x_lower,theta_y_lower,sigma_2_0_est_lower)
      est_parameters_upper=c(theta_x_upper,theta_y_upper,sigma_2_0_est_upper)

  }else if(model_name=='VFBM'){
    a_x_lower=theta[1,1]
    a_x_upper=theta[2,1]

    b_x_lower=theta[1,2]
    b_x_upper=theta[2,2]

    c_x_lower=theta[1,3]/(1+theta[1,3])
    c_x_upper=theta[2,3]/(1+theta[2,3])

    d_x_lower=theta[1,4]/(1+theta[1,4])
    d_x_upper=theta[2,4]/(1+theta[2,4])

    MSD_x_lower=rep(NA,length(d_input))
    MSD_x_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)

    #MSD_x_lower=a_x_lower*d_input^((c_x_lower*d_input)/(1+b_x_lower*d_input)+d_x_lower)
    #MSD_x_upper=a_x_upper*d_input^((c_x_upper*d_input)/(1+b_x_upper*d_input)+d_x_upper)
    if(length(index_less_than_1)>0){
      MSD_x_lower[index_less_than_1]=a_x_lower*d_input[index_less_than_1]^{(c_x_upper*d_input[index_less_than_1])/(1+b_x_upper*d_input[index_less_than_1])+d_x_upper}
      MSD_x_upper[index_less_than_1]=a_x_upper*d_input[index_less_than_1]^{(c_x_lower*d_input[index_less_than_1])/(1+b_x_lower*d_input[index_less_than_1])+d_x_lower}
      MSD_x_lower[-index_less_than_1]=a_x_lower*d_input[-index_less_than_1]^{(c_x_lower*d_input[-index_less_than_1])/(1+b_x_lower*d_input[-index_less_than_1])+d_x_lower}
      MSD_x_upper[-index_less_than_1]=a_x_upper*d_input[-index_less_than_1]^{(c_x_upper*d_input[-index_less_than_1])/(1+b_x_upper*d_input[-index_less_than_1])+d_x_upper}

    }else{
      MSD_x_lower=a_x_lower*d_input^((c_x_lower*d_input)/(1+b_x_lower*d_input)+d_x_lower)
      MSD_x_upper=a_x_upper*d_input^((c_x_upper*d_input)/(1+b_x_upper*d_input)+d_x_upper)

    }
    a_y_lower=theta[1,5]
    a_y_upper=theta[2,5]

    b_y_lower=theta[1,6]
    b_y_upper=theta[2,6]

    c_y_lower=theta[1,7]/(1+theta[1,7])
    c_y_upper=theta[2,7]/(1+theta[2,7])

    d_y_lower=theta[1,8]/(1+theta[1,8])
    d_y_upper=theta[2,8]/(1+theta[2,8])

    MSD_y_lower=rep(NA,length(d_input))
    MSD_y_upper=rep(NA,length(d_input))
    index_less_than_1=which(d_input<1)

    #MSD_y_lower=a_y_lower*d_input^((c_y_lower*d_input)/(1+b_y_lower*d_input)+d_y_lower)
    #MSD_y_upper=a_y_upper*d_input^((c_y_upper*d_input)/(1+b_y_upper*d_input)+d_y_upper)
    if(length(index_less_than_1)>0){
      MSD_y_lower[index_less_than_1]=a_y_lower*d_input[index_less_than_1]^{(c_y_upper*d_input[index_less_than_1])/(1+b_y_upper*d_input[index_less_than_1])+d_y_upper}
      MSD_y_upper[index_less_than_1]=a_y_upper*d_input[index_less_than_1]^{(c_y_lower*d_input[index_less_than_1])/(1+b_y_lower*d_input[index_less_than_1])+d_y_lower}
      MSD_y_lower[-index_less_than_1]=a_y_lower*d_input[-index_less_than_1]^{(c_y_lower*d_input[-index_less_than_1])/(1+b_y_lower*d_input[-index_less_than_1])+d_y_lower}
      MSD_y_upper[-index_less_than_1]=a_y_upper*d_input[-index_less_than_1]^{(c_y_upper*d_input[-index_less_than_1])/(1+b_y_upper*d_input[-index_less_than_1])+d_y_upper}

    }else{
      MSD_y_lower=a_y_lower*d_input^((c_y_lower*d_input)/(1+b_y_lower*d_input)+d_y_lower)
      MSD_y_upper=a_y_upper*d_input^((c_y_upper*d_input)/(1+b_y_upper*d_input)+d_y_upper)

    }
    est_parameters_lower=c(a_x_lower,b_x_lower,c_x_lower,d_x_lower,a_y_lower,
                           b_y_lower, c_y_lower,d_y_lower, sigma_2_0_est_lower);
    est_parameters_upper=c(a_x_upper,b_x_upper,c_x_upper,d_x_upper,a_y_upper,
                           b_y_upper, c_y_upper,d_y_upper, sigma_2_0_est_upper);
  }
  ans_list=list()
  ans_list$est_parameters_lower=est_parameters_lower
  ans_list$est_parameters_upper=est_parameters_upper
  ans_list$MSD_x_lower=MSD_x_lower
  ans_list$MSD_x_upper=MSD_x_upper
  ans_list$MSD_y_lower=MSD_y_lower
  ans_list$MSD_y_upper=MSD_y_upper

  return(ans_list)
}


#' Log likelihood of the model
#' @description
#' This function computes the natural logarithm of the likelihood of the
#' latent factor model for selected range of wave vectors. See 'References'.
#'
#' @param param a vector of natural logarithm of parameters
#' @param I_q_cur Fourier transformed intensity profile
#' @param B_cur current value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param index_q selected index of wave number
#' @param I_o_q_2_ori absolute square of Fourier transformed intensity profile,
#' ensemble over time
#' @param d_input sequence of lag times
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz  frame size of the intensity profile
#' @param len_t number of time steps
#' @param q wave vector in unit of um^-1
#' @param model_name model for constructing MSD, options from ('BM','OU',
#' 'FBM','OU+FBM', 'user_defined')
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return The numerical value of natural logarithm of the likelihood.
#'
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
#' @keywords internal
log_lik <- function(param,I_q_cur,B_cur,index_q,I_o_q_2_ori,d_input,
                    q_ori_ring_loc_unique_index,sz,len_t,q,model_name,
                    msd_fn=NA,msd_grad_fn=NA){

  p=length(param)-1
  theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }

  A_cur = abs(2*(I_o_q_2_ori - B_cur/2))

  ##the model is defined by MSD
  MSD_list = get_MSD_with_grad(theta,d_input,model_name, msd_fn,msd_grad_fn)
  MSD = MSD_list$msd
  log_lik_sum = 0

  NTz <- SuperGauss::NormalToeplitz$new(len_t)

  eta=B_cur/4 ##nugget

  for(i_q_selected in index_q){

    output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sqrt(sz[1]*sz[2]))
    output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sqrt(sz[1]*sz[2]))


    q_selected=q[i_q_selected]
    #beta_q = (D*q[i_q_selected]^2)
    sigma_2=A_cur[i_q_selected]/4

    acf = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
    acf[1] = acf[1]+eta
    acf=as.numeric(acf)

    log_lik_sum=log_lik_sum+sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
  }

  log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index))*log(2*pi) ##add 2pi
  if(is.nan(log_lik_sum)){
    #log_lik_sum=-10^15
    log_lik_sum=-10^50 ##make it smaller in case dealing some small value

  }
  return(log_lik_sum)
}

#' Log likelihood for anisotropic processes
#' @description
#' This function computes the natural logarithm of the likelihood of the
#' latent factor model for selected range of wave vectors of anisotropic
#' processes. See 'References'.
#'
#' @param param a vector of natural logarithm of parameters
#' @param I_q_cur Fourier transformed intensity profile
#' @param B_cur current value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param index_q selected index of wave number
#' @param I_o_q_2_ori absolute square of Fourier transformed intensity profile,
#' ensemble over time
#' @param d_input sequence of lag times
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz  frame size of the intensity profile
#' @param len_t number of time steps
#' @param q1 wave vector in unit of um^-1 in x direction
#' @param q2 wave vector in unit of um^-1 in y direction
#' @param q1_unique_index index for wave vector that give unique frequency in x direction
#' @param q2_unique_index index for wave vector that give unique frequency in y direction
#' @param model_name model for constructing MSD, options from ('BM','OU',
#' 'FBM','OU+FBM', 'user_defined')
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return The numerical value of natural logarithm of the likelihood.
#'
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
#' @keywords internal
anisotropic_log_lik <- function(param,I_q_cur,B_cur,index_q,I_o_q_2_ori,d_input,
                    q_ori_ring_loc_unique_index,sz,len_t,q1,q2,q1_unique_index,
                    q2_unique_index,model_name,msd_fn=NA,msd_grad_fn=NA){

  p=(length(param)-1)/2
  theta_x=exp(param[1:p]) ##first p parameters are parameters in ISF
  theta_y=exp(param[(p+1):(2*p)])
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[2*p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }

  A_cur=abs(2*(I_o_q_2_ori - B_cur/2))

  ##the model is defined by MSD
  MSD_list_x = get_MSD_with_grad(theta_x,d_input,model_name, msd_fn,msd_grad_fn)
  MSD_x = MSD_list_x$msd
  MSD_list_y = get_MSD_with_grad(theta_y,d_input,model_name, msd_fn,msd_grad_fn)
  MSD_y = MSD_list_y$msd

  log_lik_sum = 0

  NTz <- SuperGauss::NormalToeplitz$new(len_t)
  eta=B_cur/4 ##nugget
  q1_zero_included=c(0,q1)
  q2_zero_included=c(0,q2)

  for(i_q_selected in index_q){
    for(i_q_ori in 1:length(q_ori_ring_loc_unique_index[[i_q_selected]])){

      output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]][i_q_ori],])/(sqrt(sz[1]*sz[2]))
      output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]][i_q_ori],])/(sqrt(sz[1]*sz[2]))

      #q_selected=q[i_q_selected]
      q1_unique_index_selected=q1_unique_index[[i_q_selected]][i_q_ori]+1
      q2_unique_index_selected=q2_unique_index[[i_q_selected]][i_q_ori]+1

      sigma_2=A_cur[i_q_selected]/4

      #acf = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
      acf = sigma_2*exp(-(q1_zero_included[q1_unique_index_selected]^2*MSD_x+ q2_zero_included[q2_unique_index_selected]^2*MSD_y)/(2) )

      acf[1] = acf[1]+eta
      acf=as.numeric(acf)

      log_lik_sum=log_lik_sum+sum(NTz$logdens(z = as.numeric(output_re), acf = acf))+sum(NTz$logdens(z = as.numeric(output_im), acf = acf))
    }
  }

  log_lik_sum=log_lik_sum-0.5*sum(length(q1_unique_index)+length(q2_unique_index))*log(2*pi) ##add 2pi
  if(is.nan(log_lik_sum)){
    #log_lik_sum=-10^15
    log_lik_sum=-10^50 ##make it smaller in case dealing some small value

  }
  return(log_lik_sum)
}

#' Gradient of log likelihood
#' @description
#' This function computes the gradient for natural logarithm of the likelihood
#' for selected range of wave vectors. See 'References'.
#'
#' @param param a vector of natural logarithm of parameters
#' @param I_q_cur Fourier transformed intensity profile
#' @param B_cur current value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param index_q selected index of wave number
#' @param I_o_q_2_ori absolute square of Fourier transformed intensity profile,
#' ensemble over time
#' @param d_input sequence of lag times
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz  frame size of the intensity profile
#' @param len_t number of time steps
#' @param q wave vector in unit of um^-1
#' @param model_name stochastic process for constructing MSD, options from ('BM',
#' 'OU','FBM','OU+FBM', 'user_defined')
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return The numerical value of gradient for natural logarithm of the likelihood.
#'
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
#' @keywords internal
log_lik_grad<-function(param,I_q_cur,B_cur,index_q,I_o_q_2_ori,d_input,
                       q_ori_ring_loc_unique_index,sz,len_t,q,model_name,
                       msd_fn=NA,msd_grad_fn=NA){
  p=length(param)-1
  theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }
  #A_cur = 2*(I_o_q_2_ori - B_cur/2)
  A_cur = abs(2*(I_o_q_2_ori - B_cur/2))

  ##the model is defined by MSD
  MSD_list = get_MSD_with_grad(theta,d_input,model_name, msd_fn,msd_grad_fn)
  MSD = MSD_list$msd
  MSD_grad = MSD_list$msd_grad

  grad_trans = get_grad_trans(theta,d_input,model_name)

  eta=B_cur/4 ##nugget

  grad=rep(0,p+1)
  quad_terms=rep(0,p+1)
  trace_terms=rep(0,p+1)

  for(i_q_selected in index_q){

    output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sqrt(sz[1]*sz[2]))
    output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sqrt(sz[1]*sz[2]))

    n_q=length(q_ori_ring_loc_unique_index[[i_q_selected]])

    q_selected=q[i_q_selected]
    sigma_2=A_cur[i_q_selected]/4

    acf0 = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
    acf = acf0
    acf[1] = acf[1]+eta ##for grad this is probably no adding

    NTz=SuperGauss::Toeplitz$new(len_t, acf)
    #tilde_Sigma_inv_output_re=solve(NTz,t(output_re))
    #tilde_Sigma_inv_output_im=solve(NTz,t(output_im))

    tilde_Sigma_inv_output_re=NTz$solve(t(output_re))
    tilde_Sigma_inv_output_im=NTz$solve(t(output_im))


    acf_grad=matrix(NA,len_t,p+1)
    for(i_p in 1:p){
      acf_grad[,i_p]=-acf0*q_selected^2/4*MSD_grad[,i_p]*grad_trans[i_p]
      NTz_grad=SuperGauss::Toeplitz$new(len_t, as.numeric(acf_grad[,i_p]))
      Q_tilde_Sigma_inv_output_re=NTz_grad$prod(tilde_Sigma_inv_output_re)
      Q_tilde_Sigma_inv_output_im=NTz_grad$prod(tilde_Sigma_inv_output_im)
      quad_terms[i_p]=quad_terms[i_p]+sum(tilde_Sigma_inv_output_re*Q_tilde_Sigma_inv_output_re) ##fast way to compute quadratic terms in grad
      quad_terms[i_p]=quad_terms[i_p]+sum(tilde_Sigma_inv_output_im*Q_tilde_Sigma_inv_output_im) ##fast way to compute quadratic terms in grad

      trace_terms[i_p]=  trace_terms[i_p]+n_q*NTz$trace_grad(  as.numeric(acf_grad[,i_p]))

      #n_q*sum(diag(solve(NTz,  toeplitz(as.numeric(acf_grad[,i_p])))))
    }
    #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
    acf_grad[,p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2) ##add the sign as we use absolute value
    acf_grad[1,p+1]= acf_grad[1,p+1]+0.5
    acf_grad[,p+1]= acf_grad[,p+1]*sigma_2_0_hat ##sigma_2_0_hat is the jaccobian trans
    NTz_grad=SuperGauss::Toeplitz$new(len_t, as.numeric( acf_grad[,p+1]))
    Q_tilde_Sigma_inv_output_re=NTz_grad$prod(tilde_Sigma_inv_output_re)
    Q_tilde_Sigma_inv_output_im=NTz_grad$prod(tilde_Sigma_inv_output_im)
    quad_terms[p+1]=quad_terms[p+1]+sum(tilde_Sigma_inv_output_re*Q_tilde_Sigma_inv_output_re) ##fast way to compute quadratic terms in grad
    quad_terms[p+1]=quad_terms[p+1]+sum(tilde_Sigma_inv_output_im*Q_tilde_Sigma_inv_output_im) ##fast way to compute quadratic terms in grad

    trace_terms[p+1]=  trace_terms[p+1]+n_q*NTz$trace_grad(  as.numeric(acf_grad[,p+1]))
  }
  grad=-trace_terms+0.5*quad_terms ##note that there are two trace terms correspond to real and imaginary

  return(grad)
}

#' Gradient of log likelihood for anisotropic processes
#' @description
#' This function computes the gradient for natural logarithm of the likelihood
#' for selected range of wave vectors of anisotropic processes. See 'References'.
#'
#' @param param a vector of natural logarithm of parameters
#' @param I_q_cur Fourier transformed intensity profile
#' @param B_cur current value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param index_q selected index of wave number
#' @param I_o_q_2_ori absolute square of Fourier transformed intensity profile,
#' ensemble over time
#' @param d_input sequence of lag times
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz  frame size of the intensity profile
#' @param len_t number of time steps
#' @param q1 wave vector in unit of um^-1 in x direction
#' @param q2 wave vector in unit of um^-1 in y direction
#' @param q1_unique_index index for wave vector that give unique frequency in x direction
#' @param q2_unique_index index for wave vector that give unique frequency in y direction
#' @param model_name stochastic process for constructing MSD, options from ('BM',
#' 'OU','FBM','OU+FBM', 'user_defined')
#' @param msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return The numerical value of gradient for natural logarithm of the likelihood.
#'
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
#' @keywords internal
anisotropic_log_lik_grad<-function(param,I_q_cur,B_cur,index_q,I_o_q_2_ori,d_input,
                       q_ori_ring_loc_unique_index,sz,len_t,model_name,q1,q2,
                       q1_unique_index,q2_unique_index,msd_fn=NA,msd_grad_fn=NA){
  p=(length(param)-1)/2
  theta_x=exp(param[1:p])
  theta_y=exp(param[(p+1):(2*p)])
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[2*p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }
  #A_cur = 2*(I_o_q_2_ori - B_cur/2)
  A_cur = abs(2*(I_o_q_2_ori - B_cur/2))

  ##the model is defined by MSD
  MSD_list_x = get_MSD_with_grad(theta_x,d_input,model_name, msd_fn,msd_grad_fn)
  MSD_x = MSD_list_x$msd
  MSD_grad_x = MSD_list_x$msd_grad
  MSD_list_y = get_MSD_with_grad(theta_y,d_input,model_name, msd_fn,msd_grad_fn)
  MSD_y = MSD_list_y$msd
  MSD_grad_y = MSD_list_y$msd_grad


  grad_trans_x = get_grad_trans(theta_x,d_input,model_name)
  grad_trans_y = get_grad_trans(theta_y,d_input,model_name)

  eta=B_cur/4 ##nugget

  grad=rep(0,2*p+1)
  quad_terms=rep(0,2*p+1)
  trace_terms=rep(0,2*p+1)

  q1_zero_included=c(0,q1)
  q2_zero_included=c(0,q2)

  for(i_q_selected in index_q){
    for(i_q_ori in 1:length(q_ori_ring_loc_unique_index[[i_q_selected]])){

      output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]][i_q_ori],])/(sqrt(sz[1]*sz[2]))
      output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]][i_q_ori],])/(sqrt(sz[1]*sz[2]))

      #n_q=length(q_ori_ring_loc_unique_index[[i_q_selected]])
      n_q=1

      q1_unique_index_selected=q1_unique_index[[i_q_selected]][i_q_ori]+1
      q2_unique_index_selected=q2_unique_index[[i_q_selected]][i_q_ori]+1
      sigma_2=A_cur[i_q_selected]/4

      #acf0 = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
      acf0 = sigma_2*exp(-(q1_zero_included[q1_unique_index_selected]^2*MSD_x+
                             q2_zero_included[q2_unique_index_selected]^2*MSD_y)/(2))
      acf = acf0
      acf[1] = acf[1]+eta ##for grad this is probably no adding

      NTz=SuperGauss::Toeplitz$new(len_t, acf)
      #tilde_Sigma_inv_output_re=solve(NTz,t(output_re))
      #tilde_Sigma_inv_output_im=solve(NTz,t(output_im))

      tilde_Sigma_inv_output_re=NTz$solve(as.numeric(output_re))
      tilde_Sigma_inv_output_im=NTz$solve(as.numeric(output_im))


      acf_grad=matrix(NA,len_t,2*p+1)
      for(i_p in 1:p){
        acf_grad[,i_p]=-acf0/2*(q1_zero_included[q1_unique_index_selected]^2*MSD_grad_x[,i_p]*grad_trans_x[i_p])
        NTz_grad=SuperGauss::Toeplitz$new(len_t, as.numeric(acf_grad[,i_p]))
        Q_tilde_Sigma_inv_output_re=NTz_grad$prod(tilde_Sigma_inv_output_re)
        Q_tilde_Sigma_inv_output_im=NTz_grad$prod(tilde_Sigma_inv_output_im)
        quad_terms[i_p]=quad_terms[i_p]+sum(tilde_Sigma_inv_output_re*Q_tilde_Sigma_inv_output_re) ##fast way to compute quadratic terms in grad
        quad_terms[i_p]=quad_terms[i_p]+sum(tilde_Sigma_inv_output_im*Q_tilde_Sigma_inv_output_im) ##fast way to compute quadratic terms in grad

        trace_terms[i_p]= trace_terms[i_p]+n_q*NTz$trace_grad(as.numeric(acf_grad[,i_p]))

      #n_q*sum(diag(solve(NTz,  toeplitz(as.numeric(acf_grad[,i_p])))))
      }
      for(i_p in (p+1):(2*p)){
        acf_grad[,i_p]=-acf0/2*(q2_zero_included[q2_unique_index_selected]^2*MSD_grad_y[,(i_p-p)]*grad_trans_y[(i_p-p)])
        NTz_grad=SuperGauss::Toeplitz$new(len_t, as.numeric(acf_grad[,i_p]))
        Q_tilde_Sigma_inv_output_re=NTz_grad$prod(tilde_Sigma_inv_output_re)
        Q_tilde_Sigma_inv_output_im=NTz_grad$prod(tilde_Sigma_inv_output_im)
        quad_terms[i_p]=quad_terms[i_p]+sum(tilde_Sigma_inv_output_re*Q_tilde_Sigma_inv_output_re) ##fast way to compute quadratic terms in grad
        quad_terms[i_p]=quad_terms[i_p]+sum(tilde_Sigma_inv_output_im*Q_tilde_Sigma_inv_output_im) ##fast way to compute quadratic terms in grad

        trace_terms[i_p]= trace_terms[i_p]+n_q*NTz$trace_grad(as.numeric(acf_grad[,i_p]))

        #n_q*sum(diag(solve(NTz,  toeplitz(as.numeric(acf_grad[,i_p])))))
      }
      #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
      acf_grad[,2*p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2) ##add the sign as we use absolute value
      acf_grad[1,2*p+1]= acf_grad[1,2*p+1]+0.5
      acf_grad[,2*p+1]= acf_grad[,2*p+1]*sigma_2_0_hat ##sigma_2_0_hat is the jaccobian trans
      NTz_grad=SuperGauss::Toeplitz$new(len_t, as.numeric( acf_grad[,2*p+1]))
      Q_tilde_Sigma_inv_output_re=NTz_grad$prod(tilde_Sigma_inv_output_re)
      Q_tilde_Sigma_inv_output_im=NTz_grad$prod(tilde_Sigma_inv_output_im)
      quad_terms[2*p+1]=quad_terms[2*p+1]+sum(tilde_Sigma_inv_output_re*Q_tilde_Sigma_inv_output_re) ##fast way to compute quadratic terms in grad
      quad_terms[2*p+1]=quad_terms[2*p+1]+sum(tilde_Sigma_inv_output_im*Q_tilde_Sigma_inv_output_im) ##fast way to compute quadratic terms in grad

      trace_terms[2*p+1]= trace_terms[2*p+1]+n_q*NTz$trace_grad( as.numeric(acf_grad[,2*p+1]))
    }
  }
  grad=-trace_terms+0.5*quad_terms ##note that there are two trace terms correspond to real and imaginary

  return(grad)
}

#' Construct initial values for the parameters to be optimized over
#' @description
#' Construct initial values for the parameters to be optimized over in \code{AIUQ}
#' method of \code{SAM} class.
#'
#' @param model_name fitted model, options from  ('BM','OU','FBM','OU+FBM',
#' 'user_defined'), with Brownian motion as the default model. See 'Details'.
#' @param sigma_0_2_ini initial value for background noise, default is \code{NA}
#' @param num_param number of parameters need to be estimated in the model,
#' need to be non-NA value for \code{'user_defined'} model.
#'
#' @return A matrix with one row of initial values for the parameters to be
#' optimized over in \code{AIUQ} method of \code{SAM} class.
#' @export
#' @details
#' If \code{model_name} equals 'user_defined', then \code{num_param} need to be
#' provided to determine the length of the initial values vector.
#' @author \packageAuthor{AIUQ}
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
#' @examples
#' library(AIUQ)
#' get_initial_param(model_name = "BM")
#' @keywords internal
get_initial_param <- function(model_name,sigma_0_2_ini=NA, num_param=NA){
  if(model_name=='BM'){
    param_initial=matrix(NA,1,2) #include B
    param_initial[1,]=log(c(1,sigma_0_2_ini))#method='L-BFGS-B'
  }else if(model_name=='FBM'){
    param_initial=matrix(NA,1,3) #include B
    param_initial[1,]=log(c(rep(0.5,2),sigma_0_2_ini))#method='L-BFGS-B'
  }else if(model_name=='OU'){
    param_initial=matrix(NA,1,3) #include B
    param_initial[1,]=log(c(rep(1,2),sigma_0_2_ini))#method='L-BFGS-B'
  }else if(model_name=='OU+FBM'){
    param_initial=matrix(NA,1,5) #include B
    param_initial[1,]=log(c(rep(0.5,4),sigma_0_2_ini))#method='L-BFGS-B',
  }else if(model_name=='VFBM'){
    param_initial=matrix(NA,1,5) #include B
    param_initial[1,]=log(c(rep(0.1,4),sigma_0_2_ini))#method='L-BFGS-B}
  }else if(model_name == 'user_defined'){
    param_initial = matrix(NA,1,(num_param+1)) #include B
    param_initial[1,] = log(c(rep(0.5,num_param),sigma_0_2_ini))
  }else if(model_name=='VFBM_anisotropic'){
      param_initial=matrix(NA,1,9) #include B
      param_initial[1,]=log(c(rep(0.1,8),sigma_0_2_ini))#method='L-BFGS-B'
  }else if(model_name=='BM_anisotropic'){
    param_initial=matrix(NA,1,3) #include B
    param_initial[1,]=log(c(1,1,sigma_0_2_ini))#method='L-BFGS-B',
  }else if(model_name=='FBM_anisotropic'){
    param_initial=matrix(NA,1,5) #include B
    param_initial[1,]=log(c(0.5,2,0.5,2,sigma_0_2_ini))#method='L-BFGS-B',
  }else if(model_name=='OU_anisotropic'){
    param_initial=matrix(NA,1,5) #include B
    param_initial[1,]=log(c(rep(1,4),sigma_0_2_ini))#method='L-BFGS-B',
  }else if(model_name=='OU+BM_anisotropic'){
    param_initial=matrix(NA,1,7) #include B
    param_initial[1,]=log(c(rep(0.5,6),sigma_0_2_ini))#method='L-BFGS-B',
  }else if(model_name=='OU+FBM_anisotropic'){
    param_initial=matrix(NA,1,9) #include B
    param_initial[1,]=log(c(rep(c(1,1,0.5,2),2),sigma_0_2_ini))#method='L-BFGS-B',
  }else if(model_name == 'user_defined_anisotropic'){
    param_initial = matrix(NA,1,2*num_param+1) #include B
    param_initial[1,] = log(c(rep(0.5,2*num_param),sigma_0_2_ini))
  }
  return(param_initial)
}

#' Construct parameter transformation for optimization using exact gradient
#' @description
#' Construct parameter transformation for parameters to be optimized over in \code{AIUQ}
#' method of \code{SAM} class. See 'References'.
#'
#' @param theta parameters to be optimized over
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#'
#' @return A vector of transformed parameters to be optimized over in \code{AIUQ}
#' method of \code{SAM} class.
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
#' @keywords internal
get_grad_trans<-function(theta,d_input,model_name){
  if(model_name=='BM'||model_name=='BM_anisotropic'){
    grad_trans = theta[1]
  }else if(model_name=='FBM'||model_name=='FBM_anisotropic'){
    beta = theta[1]
    alpha = 2*theta[2]/(1+theta[2])

    grad_trans = c(beta, alpha*(1-alpha/2))
  }else if(model_name=='OU'||model_name=='OU_anisotropic'){
    rho = theta[1]/(1+theta[1])
    amplitude = theta[2]

    grad_trans = c(rho*(1-rho), amplitude)
  }else if(model_name=='OU+FBM'||model_name=='OU+FBM_anisotropic'){
    rho = theta[1]/(1+theta[1])
    amplitude = theta[2]
    beta = theta[3]
    alpha = 2*theta[4]/(1+theta[4])

    grad_trans = c(rho*(1-rho), amplitude, beta, alpha*(1-alpha/2))
  }else if(model_name=='user_defined'||model_name=='user_defined_anisotropic'){
    grad_trans = theta
  }else if(model_name=='VFBM'||model_name=='VFBM_anisotropic'){
    a = theta[1]
    b = theta[2]
    c = theta[3]/(1+theta[3])
    d = theta[4]/(1+theta[4])

    #grad_trans = c(a,b,c*(1-c),d*(1-d))
    grad_trans = c(a,b,c/(1-c),d/(1-d))
  }
  return(grad_trans)
}

#' Compute 95% confidence interval
#' @description
#' This function construct the lower and upper bound for 95% confidence interval
#' of estimated parameters for the given model, including parameters contained
#' in the intermediate scattering function and background noise. See 'References'.
#'
#' @param param_est a vector of natural logarithm of estimated parameters from
#' maximize the log likelihood. This vector will serve as initial values in the
#' \code{optim} function.
#' @param I_q_cur Fourier transformed intensity profile
#' @param B_cur current value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param index_q selected index of wave number
#' @param I_o_q_2_ori absolute square of Fourier transformed intensity profile,
#' ensemble over time
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz  frame size of the intensity profile
#' @param len_t number of time steps
#' @param q wave vector in unit of um^-1
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#' @param estimation_method method for constructing 95% confidence interval,
#' default is asymptotic
#' @param M number of particles
#' @param num_iteration_max the maximum number of iterations in \code{optim}
#' @param lower_bound lower bound for the "L-BFGS-B" method in \code{optim}
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return A matrix of lower and upper bound for natural logarithm of
#' parameters in the fitted model using \code{AIUQ} method in \code{SAM} class
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
#' @keywords internal
param_uncertainty<-function(param_est,I_q_cur,B_cur=NA,index_q,
                            I_o_q_2_ori,q_ori_ring_loc_unique_index,
                            sz,len_t,d_input,q,model_name,
                            estimation_method='asymptotic',M,
                            num_iteration_max,lower_bound,msd_fn=NA,
                            msd_grad_fn=NA){
  p = length(param_est)-1
  q_lower=q-min(q)
  if(model_name == "user_defined"){
    if(is.function(msd_grad_fn)==T){
      gr = log_lik_grad
    }else{gr = NULL}
  }else{gr = log_lik_grad}
  #param_ini=param_est
  m_param_lower = try(optim(param_est,log_lik,gr = gr,I_q_cur=I_q_cur,B_cur=NA,
                            index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                            q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                            method='L-BFGS-B',lower=lower_bound,
                            control = list(fnscale=-1,maxit=num_iteration_max),
                            sz=sz,len_t=len_t,d_input=d_input,q=q_lower,
                            model_name=model_name,msd_fn=msd_fn,
                            msd_grad_fn=msd_grad_fn),TRUE)


  if(class(m_param_lower)[1]=="try-error"){
    compute_twice=T
    m_param_lower = try(optim(param_est+c(rep(0.5,p),0),log_lik,gr=gr,
                              I_q_cur=I_q_cur,B_cur=NA,
                              index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                              q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                              method='L-BFGS-B',lower=lower_bound,
                              control = list(fnscale=-1,maxit=num_iteration_max),
                              sz=sz,len_t=len_t,d_input=d_input,q=q_lower,
                              model_name=model_name,msd_fn=msd_fn,
                              msd_grad_fn=msd_grad_fn),TRUE)

  }
  q_upper=q+min(q)

  m_param_upper = try(optim(param_est,log_lik,gr=gr,
                            I_q_cur=I_q_cur,B_cur=NA,
                            index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                            q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                            method='L-BFGS-B',lower=lower_bound,
                            control = list(fnscale=-1,maxit=num_iteration_max),
                            sz=sz,len_t=len_t,d_input=d_input,q=q_upper,
                            model_name=model_name,msd_fn=msd_fn,
                            msd_grad_fn=msd_grad_fn),TRUE)

  if(class(m_param_upper)[1]=="try-error"){
    compute_twice = T
    m_param_upper = try(optim(param_est+c(rep(0.5,p),0),log_lik,gr=gr,
                              I_q_cur=I_q_cur,B_cur=NA,
                              index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                              q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                              method='L-BFGS-B',lower=lower_bound,
                              control = list(fnscale=-1,maxit=num_iteration_max),
                              sz=sz,len_t=len_t,d_input=d_input,q=q_upper,
                              model_name=model_name,msd_fn=msd_fn,
                              msd_grad_fn=msd_grad_fn),TRUE)
  }

  param_range = matrix(NA,2,p+1)
  for(i in 1:(p+1) ){
    param_range[1,i] = min(m_param_lower$par[i],m_param_upper$par[i])
    param_range[2,i] = max(m_param_lower$par[i],m_param_upper$par[i])

  }
  half_length_param_range_fft = (param_range[2,]-param_range[1,])/2

  if(estimation_method=='asymptotic'){
    theta = exp(param_est[-(p+1)]) ##first p parameters are parameters in ISF
    if(is.na(B_cur)){ ##this fix the dimension
      sigma_2_0_hat = exp(param_est[p+1]) ##noise
      B_cur = 2*sigma_2_0_hat
    }

    A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
    eta = B_cur/4 ##nugget

    MSD_list = get_MSD_with_grad(theta,d_input,model_name,msd_fn,msd_grad_fn)
    MSD = MSD_list$msd
    MSD_grad = MSD_list$msd_grad

    grad_trans = get_grad_trans(theta,d_input,model_name)
    Hessian_list = as.list(index_q)
    Hessian_sum = 0

    for(i_q_selected in index_q){
      q_selected=q[i_q_selected]

      sigma_2=A_cur[i_q_selected]/4

      acf0=sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
      acf=acf0
      acf[1] = acf[1]+eta ##for grad this is probably no adding
      acf=as.numeric(acf)

      Tz <- SuperGauss::Toeplitz$new(len_t,acf=acf)
      Hessian = matrix(NA,p+1,p+1) ##last one is
      acf_grad = matrix(NA,len_t,p+1)
      for(i_p in 1:p){
        acf_grad[,i_p] = -acf0*q_selected^2/4* MSD_grad[,i_p]*grad_trans[i_p]
      }
      #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
      acf_grad[,p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2)
      acf_grad[1,p+1]= acf_grad[1,p+1]+0.5
      acf_grad[,p+1]= acf_grad[,p+1]*sigma_2_0_hat

      for(i_p in 1:(p+1) ){
        for(j_p in 1:(p+1) ){
          Hessian[i_p,j_p]=Tz$trace_hess(as.numeric(acf_grad[,i_p]), as.numeric(acf_grad[,j_p]) )
        }
      }

      #Hessian_list[[i_q_selected]]=Hessian
      Hessian_sum=Hessian_sum+Hessian*length(q_ori_ring_loc_unique_index[[i_q_selected]])
      ###a litte more conservation is to say they are perfectly correlated in a ring
      #Hessian_sum=Hessian_sum+Hessian
    }

    Hessian_sum = Hessian_sum*M/sum(lengths(q_ori_ring_loc_unique_index[index_q]))
    if(kappa(Hessian_sum)>1e10){
      epsilon <- 1e-6
      Hessian_sum <- Hessian_sum + epsilon * diag(ncol(Hessian_sum))
    }
    sd_theta_B = sqrt(diag(solve(Hessian_sum)))
    #sd_theta_B=sqrt(diag(Hessian_inv_sum/sum(lengths(q_ori_ring_loc_unique_index))^2 ))
    param_range[1,]=param_range[1,]-sd_theta_B*qnorm(0.975) ##this is 10 times larger to account for not estimating A and B correct and other misspecification
    param_range[2,]=param_range[2,]+sd_theta_B*qnorm(0.975)
    half_length_param_range_est = sd_theta_B*qnorm(0.975)
  }
  return(param_range)
}

#' Compute 95% confidence interval for anisotropic processes
#' @description
#' This function construct the lower and upper bound for 95% confidence interval
#' of estimated parameters for the given anisotropic model, including parameters
#' contained in the intermediate scattering function and background noise.
#' See 'References'.
#'
#' @param param_est a vector of natural logarithm of estimated parameters from
#' maximize the log likelihood. This vector will serve as initial values in the
#' \code{optim} function.
#' @param I_q_cur Fourier transformed intensity profile
#' @param B_cur current value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param index_q selected index of wave number
#' @param I_o_q_2_ori absolute square of Fourier transformed intensity profile,
#' ensemble over time
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz  frame size of the intensity profile
#' @param len_t number of time steps
#' @param q1 wave vector in unit of um^-1 in x direction
#' @param q2 wave vector in unit of um^-1 in y direction
#' @param q1_unique_index index for wave vector that give unique frequency in x direction
#' @param q2_unique_index index for wave vector that give unique frequency in y direction
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#' @param estimation_method method for constructing 95% confidence interval,
#' default is asymptotic
#' @param M number of particles
#' @param num_iteration_max the maximum number of iterations in \code{optim}
#' @param lower_bound lower bound for the "L-BFGS-B" method in \code{optim}
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return A matrix of lower and upper bound for natural logarithm of
#' parameters in the fitted model using \code{AIUQ} method in \code{SAM} class
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
#' @keywords internal
param_uncertainty_anisotropic<-function(param_est,I_q_cur,B_cur=NA,index_q,
                            I_o_q_2_ori,q_ori_ring_loc_unique_index,
                            sz,len_t,d_input,q1,q2,q1_unique_index,q2_unique_index,
                            model_name,estimation_method='asymptotic',M,
                            num_iteration_max,lower_bound,msd_fn=NA,
                            msd_grad_fn=NA){
  p=(length(param_est)-1)/2
  q1_lower=q1-min(q1)
  q2_lower=q2-min(q2)

  # if(model_name == "user_defined"){
  #   if(is.function(msd_grad_fn)==T){
  #     gr = log_lik_grad
  #   }else{gr = NULL}
  # }else{gr = log_lik_grad}

  #param_ini=param_est
  m_param_lower = try(optim(param_est,anisotropic_log_lik,#gr = gr,
                            I_q_cur=I_q_cur,B_cur=NA,
                            index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                            q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                            method='L-BFGS-B',lower=lower_bound,
                            control = list(fnscale=-1,maxit=num_iteration_max),
                            sz=sz,len_t=len_t,d_input=d_input,q1=q1_lower,
                            q2=q2_lower, q1_unique_index=q1_unique_index,
                            q2_unique_index=q2_unique_index,
                            model_name=model_name,msd_fn=msd_fn,
                            msd_grad_fn=msd_grad_fn),TRUE)


  if(class(m_param_lower)[1]=="try-error"){
    compute_twice=T
    m_param_lower = try(optim(param_est,anisotropic_log_lik,#gr = gr,
                              I_q_cur=I_q_cur,B_cur=NA,
                              index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                              q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                              method='L-BFGS-B',lower=lower_bound,
                              control = list(fnscale=-1,maxit=num_iteration_max),
                              sz=sz,len_t=len_t,d_input=d_input,q1=q1_lower,
                              q2=q2_lower, q1_unique_index=q1_unique_index,
                              q2_unique_index=q2_unique_index,
                              model_name=model_name,msd_fn=msd_fn,
                              msd_grad_fn=msd_grad_fn),TRUE)

  }
  q1_upper=q1+min(q1)
  q2_upper=q2+min(q2)

  m_param_upper = try(optim(param_est,anisotropic_log_lik, #gr=gr,
                            I_q_cur=I_q_cur,B_cur=NA,
                            index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                            q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                            method='L-BFGS-B',lower=lower_bound,
                            control = list(fnscale=-1,maxit=num_iteration_max),
                            sz=sz,len_t=len_t,d_input=d_input,q1=q1_upper,
                            q2=q2_upper, q1_unique_index=q1_unique_index,
                            q2_unique_index=q2_unique_index,
                            model_name=model_name,msd_fn=msd_fn,
                            msd_grad_fn=msd_grad_fn),TRUE)

  if(class(m_param_upper)[1]=="try-error"){
    compute_twice = T
    m_param_upper = try(optim(param_est+c(rep(0.5,p),0),anisotropic_log_lik,#gr=gr,
                              I_q_cur=I_q_cur,B_cur=NA,
                              index_q=index_q,I_o_q_2_ori=I_o_q_2_ori,
                              q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                              method='L-BFGS-B',lower=lower_bound,
                              control = list(fnscale=-1,maxit=num_iteration_max),
                              sz=sz,len_t=len_t,d_input=d_input,q1=q1_upper,
                              q2=q2_upper, q1_unique_index=q1_unique_index,
                              q2_unique_index=q2_unique_index,
                              model_name=model_name,msd_fn=msd_fn,
                              msd_grad_fn=msd_grad_fn),TRUE)
  }

  param_range = matrix(NA,2,2*p+1)
  for(i in 1:(2*p+1) ){
    param_range[1,i] = min(m_param_lower$par[i],m_param_upper$par[i])
    param_range[2,i] = max(m_param_lower$par[i],m_param_upper$par[i])

  }
  half_length_param_range_fft = (param_range[2,]-param_range[1,])/2

  if(estimation_method=='asymptotic'){
    theta_x = exp(param_est[1:p]) ##first p parameters are parameters in ISF
    theta_y = exp(param_est[(p+1):(2*p)])
    #theta = exp(param_est[-(p+1)])
    if(is.na(B_cur)){ ##this fix the dimension
      sigma_2_0_hat = exp(param_est[2*p+1]) ##noise
      B_cur = 2*sigma_2_0_hat
    }

    A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
    eta = B_cur/4 ##nugget

    MSD_list_x = get_MSD_with_grad(theta_x,d_input,model_name,msd_fn,msd_grad_fn)
    MSD_x = MSD_list_x$msd
    MSD_grad_x = MSD_list_x$msd_grad
    MSD_list_y = get_MSD_with_grad(theta_y,d_input,model_name,msd_fn,msd_grad_fn)
    MSD_y = MSD_list_y$msd
    MSD_grad_y = MSD_list_y$msd_grad

    grad_trans_x = get_grad_trans(theta_x,d_input,model_name)
    grad_trans_y = get_grad_trans(theta_y,d_input,model_name)
    Hessian_list = as.list(index_q)
    Hessian_sum = 0

    q1_zero_included=c(0,q1)
    q2_zero_included=c(0,q2)

    for(i_q_selected in index_q){
      for(i_q_ori in 1:length(q_ori_ring_loc_unique_index[[i_q_selected]]) ){
        #q_selected=q[i_q_selected]
        q1_unique_index_selected=q1_unique_index[[i_q_selected]][i_q_ori]+1
        q2_unique_index_selected=q2_unique_index[[i_q_selected]][i_q_ori]+1

        sigma_2=A_cur[i_q_selected]/4

        #acf0=sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
        acf0 = sigma_2*exp(-(q1_zero_included[q1_unique_index_selected]^2*MSD_x+
                               q2_zero_included[q2_unique_index_selected]^2*MSD_y)/(2) ) ##assume 2d
        acf=acf0
        acf[1] = acf[1]+eta ##for grad this is probably no adding
        acf=as.numeric(acf)

        Tz <- SuperGauss::Toeplitz$new(len_t,acf=acf)
        Hessian = matrix(NA,2*p+1,2*p+1) ##last one is
        acf_grad = matrix(NA,len_t,2*p+1)
        for(i_p in 1:p){
          #acf_grad[,i_p] = -acf0*q_selected^2/4* MSD_grad[,i_p]*grad_trans[i_p]
          acf_grad[,i_p]=-acf0*q1_zero_included[q1_unique_index_selected]^2/2*
            MSD_grad_x[,i_p]*grad_trans_x[i_p]
          acf_grad[,p+i_p]=-acf0*q2_zero_included[q2_unique_index_selected]^2/2*
            MSD_grad_y[,i_p]*grad_trans_y[i_p]
        }
        #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
        acf_grad[,2*p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2)
        acf_grad[1,2*p+1]= acf_grad[1,2*p+1]+0.5
        acf_grad[,2*p+1]= acf_grad[,2*p+1]*sigma_2_0_hat

        for(i_p in 1:(2*p+1) ){
          for(j_p in 1:(2*p+1) ){
            Hessian[i_p,j_p]=Tz$trace_hess(as.numeric(acf_grad[,i_p]), as.numeric(acf_grad[,j_p]) )
          }
        }

        #Hessian_list[[i_q_selected]]=Hessian
        Hessian_sum=Hessian_sum+Hessian #length(q_ori_ring_loc_unique_index[[i_q_selected]])
        ###a litte more conservation is to say they are perfectly correlated in a ring
        #Hessian_sum=Hessian_sum+Hessian
      }
    }

    Hessian_sum = Hessian_sum*M/sum(lengths(q_ori_ring_loc_unique_index[index_q]))
    if(kappa(Hessian_sum)>1e10){
      epsilon <- 1e-6
      Hessian_sum <- Hessian_sum + epsilon * diag(ncol(Hessian_sum))
    }
    sd_theta_B = sqrt(diag(solve(Hessian_sum)))
    #sd_theta_B=sqrt(diag(Hessian_inv_sum/sum(lengths(q_ori_ring_loc_unique_index))^2 ))
    param_range[1,]=param_range[1,]-sd_theta_B*qnorm(0.975) ##this is 10 times larger to account for not estimating A and B correct and other misspecification
    param_range[2,]=param_range[2,]+sd_theta_B*qnorm(0.975)
    half_length_param_range_est = sd_theta_B*qnorm(0.975)
  }
  return(param_range)
}

#' Simulate 2D particle trajectory follows Brownian Motion
#'
#' @description
#' Simulate 2D particle trajectory follows Brownian Motion (BM) for \code{M}
#' particles.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma distance moved per time step
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma = 0.5
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = bm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma)
#' @keywords internal
bm_particle_intensity <- function(pos0,M,len_t,sigma){
  pos = matrix(NA,M*len_t,2)
  pos[1:M,] = pos0
  for(i in 1:(len_t-1)){
    pos[i*M+(1:M),] = pos[(i-1)*M+(1:M),]+matrix(rnorm(2*M,sd=sigma),M,2)
  }
  return(pos)
}

#' Simulate 2D particle trajectory follows anisotropic Brownian Motion
#'
#' @description
#' Simulate 2D particle trajectory follows anisotropic Brownian Motion (BM) for
#' \code{M} particles, with different step sizes in x, y-directions.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma distance moved per time step in x,y-directions, a vector of length 2
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma = c(0.5,0.1)
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = anisotropic_bm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma)
#' @keywords internal
anisotropic_bm_particle_intensity <- function(pos0,M,len_t,sigma){
  pos = matrix(NA,M*len_t,2)
  pos[1:M,] = pos0
  for(i in 1:(len_t-1)){
    pos[i*M+(1:M),1] = pos[(i-1)*M+(1:M),1]+matrix(rnorm(M,sd=sigma[1]),M,1)
    pos[i*M+(1:M),2] = pos[(i-1)*M+(1:M),2]+matrix(rnorm(M,sd=sigma[2]),M,1)
  }
  return(pos)
}


#' Simulate 2D particle trajectory follows OU process
#' @description
#' Simulate 2D particle trajectory follows Ornstein–Uhlenbeck process(OU) for
#' \code{M} particles.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma distance moved per time step
#' @param rho correlation between successive step and previous step,
#' value between 0 and 1
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma = 2
#' rho = 0.95
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = ou_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma, rho=rho)
#' @keywords internal
ou_particle_intensity <- function(pos0,M,len_t,sigma,rho){
  pos = matrix(NA,M*len_t,2)
  pos[1:M,] = pos0+matrix(rnorm(2*M, sd=sigma),M,2)
  sd_innovation_OU = sqrt(sigma^2*(1-rho^2))

  for(i in 1:(len_t-1)){
    pos[i*M+(1:M),1] = rho*(pos[(i-1)*M+(1:M),1]-pos0[,1])+pos0[,1]+sd_innovation_OU*rnorm(M)
    pos[i*M+(1:M),2] = rho*(pos[(i-1)*M+(1:M),2]-pos0[,2])+pos0[,2]+sd_innovation_OU*rnorm(M)
  }
  return(pos)
}

#' Simulate 2D particle trajectory follows anisotropic OU process
#' @description
#' Simulate 2D particle trajectory follows anisotropic Ornstein–Uhlenbeck
#' process(OU) for \code{M} particles, with different step sizes in x, y-directions.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma distance moved per time step in x, y-directions, a vector of length 2
#' @param rho correlation between successive step and previous step in x, y-directions,
#' a vector of length 2 with values between 0 and 1
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma = c(2,2.5)
#' rho = c(0.95,0.9)
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = anisotropic_ou_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma, rho=rho)
#' @keywords internal
anisotropic_ou_particle_intensity <- function(pos0,M,len_t,sigma,rho){
  pos = matrix(NA,M*len_t,2)
  pos[1:M,1] = pos0[,1]+matrix(rnorm(M, sd=sigma[1]),M,1)
  pos[1:M,2] = pos0[,2]+matrix(rnorm(M, sd=sigma[2]),M,1)
  sd_innovation_OU_1 = sqrt(sigma[1]^2*(1-rho[1]^2))
  sd_innovation_OU_2 = sqrt(sigma[2]^2*(1-rho[2]^2))

  for(i in 1:(len_t-1)){
    pos[i*M+(1:M),1] = rho[1]*(pos[(i-1)*M+(1:M),1]-pos0[,1])+pos0[,1]+sd_innovation_OU_1*matrix(rnorm(M),M,1)
    pos[i*M+(1:M),2] = rho[2]*(pos[(i-1)*M+(1:M),2]-pos0[,2])+pos0[,2]+sd_innovation_OU_2*matrix(rnorm(M),M,1)
  }
  return(pos)
}

#' Construct correlation matrix for FBM
#' @description
#' Construct correlation matrix for fractional Brownian motion.
#'
#' @param len_t number of time steps
#' @param H Hurst parameter, value between 0 and 1
#'
#' @return Correlation matrix with dimension \code{len_t-1} by \code{len_t-1}.
#' @export
#' @author \packageAuthor{AIUQ}
#'
#' @examples
#' library(AIUQ)
#' len_t = 50
#' H = 0.3
#' m = corr_fBM(len_t=len_t,H=H)
#' @keywords internal
# corr_fBM <- function(len_t,H){
#   corr = matrix(NA, len_t, len_t)
#   for(i in 1:(len_t)){
#     for(j in 1:(len_t)){
#       corr[i,j] = 0.5*(i^(2*H)+j^(2*H)-abs(i-j)^(2*H))
#     }
#   }
#   return(corr)
# }
corr_fBM <- function(len_t,H){
  cov = matrix(NA, len_t-1, len_t-1)
  for(i in 0:(len_t-2)){
  cov[i+1,] = 0.5*(abs(i-(0:(len_t-2))+1)^(2*H))+0.5*(abs(1-i+(0:(len_t-2)))^(2*H))-abs(i-(0:(len_t-2)))^(2*H)
  }
  return(cov)
}

#' Simulate 2D particle trajectory follows fBM
#' @description
#' Simulate 2D particle trajectory follows fraction Brownian Motion(fBM) for
#' \code{M} particles.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma distance moved per time step
#' @param H Hurst parameter, value between 0 and 1
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma = 2
#' H = 0.3
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = fbm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma,H=H)
#' @keywords internal
# fbm_particle_intensity <- function(pos0,M,len_t,sigma,H){
#   pos = matrix(NA,M*len_t,2)
#   pos[,1] = rep(pos0[,1],len_t)
#   pos[,2] = rep(pos0[,2],len_t)
#   fBM_corr = corr_fBM(len_t,H)
#   L = t(chol(sigma^2*fBM_corr))
#   pos[,1] = pos[,1]+as.numeric(t(L%*%matrix(rnorm((len_t)*M),nrow=len_t,ncol=M)))
#   pos[,2] = pos[,2]+as.numeric(t(L%*%matrix(rnorm((len_t)*M),nrow=len_t,ncol=M)))
#   return(pos)
# }
fbm_particle_intensity <- function(pos0,M,len_t,sigma,H){
  pos = matrix(NA,M*len_t,2)
  pos[,1] = rep(pos0[,1],len_t)
  pos[,2] = rep(pos0[,2],len_t)
  fBM_corr = corr_fBM(len_t,H)
  L = t(chol(sigma^2*fBM_corr))
  increments1 = L%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  increments2 = L%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  pos[(M+1):(M*len_t),1] = pos[(M+1):(M*len_t),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos[(M+1):(M*len_t),2] = pos[(M+1):(M*len_t),2]+as.numeric(t(apply(increments2,2,cumsum)))
  return(pos)
}

#' Simulate 2D particle trajectory follows anisotropic fBM
#' @description
#' Simulate 2D particle trajectory follows anisotropic fraction Brownian Motion(fBM) for
#' \code{M} particles, with different step sizes in x, y-directions.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma distance moved per time step in x, y-directions, a vector of length 2
#' @param H Hurst parameter in x, y-directions, a vector of length 2, value
#' between 0 and 1
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma = c(2,1)
#' H = c(0.3,0.4)
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = anisotropic_fbm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma,H=H)
#' @keywords internal
anisotropic_fbm_particle_intensity <- function(pos0,M,len_t,sigma,H){
  pos = matrix(NA,M*len_t,2)
  pos[,1] = rep(pos0[,1],len_t)
  pos[,2] = rep(pos0[,2],len_t)
  fBM_corr1 = corr_fBM(len_t,H[1])
  L1 = t(chol(sigma[1]^2*fBM_corr1))
  fBM_corr2 = corr_fBM(len_t,H[2])
  L2 = t(chol(sigma[2]^2*fBM_corr2))
  increments1 = L1%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  increments2 = L2%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  pos[(M+1):(M*len_t),1] = pos[(M+1):(M*len_t),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos[(M+1):(M*len_t),2] = pos[(M+1):(M*len_t),2]+as.numeric(t(apply(increments2,2,cumsum)))
  return(pos)
}

#' Simulate 2D particle trajectory follows fBM plus OU
#' @description
#' Simulate 2D particle trajectory follows fraction Brownian Motion(fBM) plus a
#' Ornstein–Uhlenbeck(OU) process for \code{M} particles.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma_fbm distance moved per time step in fractional Brownian Motion
#' @param sigma_ou distance moved per time step in Ornstein–Uhlenbeck process
#' @param H Hurst parameter of fractional Brownian Motion, value between 0 and 1
#' @param rho correlation between successive step and previous step in OU process,
#' value between 0 and 1
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma_fbm = 2
#' H = 0.3
#' sigma_ou = 2
#' rho = 0.95
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#'
#' pos = fbm_ou_particle_intensity(pos0=pos0, M=M, len_t=len_t,
#'   sigma_fbm=sigma_fbm, sigma_ou=sigma_ou, H=H, rho=rho)
#' @keywords internal
# fbm_ou_particle_intensity <- function(pos0,M,len_t,sigma_fbm,sigma_ou,H,rho){
#   pos_ou = matrix(NA,M*len_t,2)
#   pos_fbm = matrix(NA,M*len_t,2)
#
#   pos0 = pos0 + matrix(rnorm(2*M, sd=sigma_ou),M,2)
#   pos_ou[1:M,] = pos0
#   sd_innovation_OU=sqrt(sigma_ou^2*(1-rho^2))
#
#   for(i in 1:(len_t-1)){
#     pos_ou[i*M+(1:M),] = rho*(pos_ou[(i-1)*M+(1:M),]-pos0)+pos0+
#       sd_innovation_OU*matrix(rnorm(2*M),M,2)
#   }
#
#
#   fBM_corr = corr_fBM(len_t,H)
#   L = t(chol(sigma_fbm^2*fBM_corr))
#   pos_fbm[,1] = as.numeric(t(L%*%matrix(rnorm((len_t)*M),nrow=len_t,ncol=M)))
#   pos_fbm[,2] = as.numeric(t(L%*%matrix(rnorm((len_t)*M),nrow=len_t,ncol=M)))
#   pos = pos_ou+pos_fbm
#
#   return(pos)
# }
fbm_ou_particle_intensity <- function(pos0,M,len_t,sigma_fbm,sigma_ou,H,rho){
  pos1 = matrix(NA,M*len_t,2)
  pos2 = matrix(NA,M*len_t,2)
  pos = matrix(NA,M*len_t,2)
  pos0 = pos0 + matrix(rnorm(2*M, sd=sigma_ou),M,2)
  pos[1:M,] = pos0
  pos1[1:M,] = pos0
  sd_innovation_OU=sqrt(sigma_ou^2*(1-rho^2))

  for(i in 1:(len_t-1)){
    pos1[i*M+(1:M),] = rho*(pos1[(i-1)*M+(1:M),]-pos0)+pos0+sd_innovation_OU*matrix(rnorm(2*M),M,2)
  }

  pos2[,1] = rep(pos0[,1],len_t)
  pos2[,2] = rep(pos0[,2],len_t)
  fBM_corr = corr_fBM(len_t,H)
  L = t(chol(sigma_fbm^2*fBM_corr))
  increments1 = L%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  increments2 = L%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  pos2[(M+1):(M*len_t),1] = pos2[(M+1):(M*len_t),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos2[(M+1):(M*len_t),2] = pos2[(M+1):(M*len_t),2]+as.numeric(t(apply(increments2,2,cumsum)))
  pos = pos1+pos2 - cbind(rep(pos0[,1],len_t), rep(pos0[,2],len_t))

  return(pos)
}


#' Simulate 2D particle trajectory follows anisotropic fBM plus OU
#' @description
#' Simulate 2D particle trajectory follows anisotropic fraction Brownian
#' Motion(fBM) plus a Ornstein–Uhlenbeck(OU) process for \code{M} particles,
#' with different step sizes in x, y-directions.
#'
#' @param pos0 initial position for \code{M} particles, matrix with dimension M by 2
#' @param M number of particles
#' @param len_t number of time steps
#' @param sigma_fbm distance moved per time step in fractional Brownian Motion
#' in x, y-directions, a vector of length 2
#' @param sigma_ou distance moved per time step in Ornstein–Uhlenbeck process
#' in x, y-directions, a vector of length 2
#' @param H Hurst parameter of fractional Brownian Motion in x, y-directions,
#' a vector of length 2 with values between 0 and 1
#' @param rho correlation between successive step and previous step in OU process
#' in x, y-directions, a vector of length 2 with values between 0 and 1
#'
#' @return Position matrix with dimension \code{M}\eqn{\times}{%\times}\code{len_t}
#' by 2 for particle trajectory. The first \code{M} rows being the initial position
#' \code{pos0}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' M = 10
#' len_t = 50
#' sigma_fbm = c(2,1)
#' H = c(0.3,0.4)
#' sigma_ou = c(2,2.5)
#' rho = c(0.95,0.9)
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#'
#' pos = anisotropic_fbm_ou_particle_intensity(pos0=pos0, M=M, len_t=len_t,
#'   sigma_fbm=sigma_fbm, sigma_ou=sigma_ou, H=H, rho=rho)
#' @keywords internal
anisotropic_fbm_ou_particle_intensity <- function(pos0,M,len_t,sigma_fbm,sigma_ou,H,rho){
  pos1 = matrix(NA,M*len_t,2)
  pos2 = matrix(NA,M*len_t,2)
  pos = matrix(NA,M*len_t,2)
  pos0[,1] = pos0[,1] + matrix(rnorm(M, sd=sigma_ou[1]),M,1)
  pos0[,2] = pos0[,2] + matrix(rnorm(M, sd=sigma_ou[2]),M,1)

  pos[1:M,] = pos0
  pos1[1:M,] = pos0
  sd_innovation_OU_1=sqrt(sigma_ou[1]^2*(1-rho[1]^2))
  sd_innovation_OU_2=sqrt(sigma_ou[2]^2*(1-rho[2]^2))

  for(i in 1:(len_t-1)){
    pos1[i*M+(1:M),1] = rho*(pos1[(i-1)*M+(1:M),1]-pos0[,1])+pos0[,1]+sd_innovation_OU_1*matrix(rnorm(M),M,1)
    pos1[i*M+(1:M),2] = rho*(pos1[(i-1)*M+(1:M),2]-pos0[,2])+pos0[,2]+sd_innovation_OU_2*matrix(rnorm(M),M,1)
  }

  pos2[,1] = rep(pos0[,1],len_t)
  pos2[,2] = rep(pos0[,2],len_t)

  fBM_corr1 = corr_fBM(len_t,H[1])
  L1 = t(chol(sigma_fbm[1]^2*fBM_corr1))
  fBM_corr2 = corr_fBM(len_t,H[2])
  L2 = t(chol(sigma_fbm[2]^2*fBM_corr2))

  increments1 = L1%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  increments2 = L2%*%matrix(rnorm((len_t-1)*M),nrow=len_t-1,ncol=M)
  pos2[(M+1):(M*len_t),1] = pos2[(M+1):(M*len_t),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos2[(M+1):(M*len_t),2] = pos2[(M+1):(M*len_t),2]+as.numeric(t(apply(increments2,2,cumsum)))
  pos = pos1+pos2 - cbind(rep(pos0[,1],len_t), rep(pos0[,2],len_t))

  return(pos)
}


#' Construct intensity profile for a given particle trajectory
#' @description
#' Construct intensity profile with structure 'T_SS_mat' for a given particle
#' trajectory, background intensity profile, and user defined radius of particle.
#'
#' @param len_t number of time steps
#' @param M number of particles
#' @param I background intensity profile. See 'Details'.
#' @param pos position matrix for particle trajectory
#' @param Ic vector of maximum intensity of each particle
#' @param sz frame size of simulated square image
#' @param sigma_p radius of the spherical particle (3sigma_p)
#'
#' @return Intensity profile matrix with structure 'T_SS_mat' (matrix with
#' dimension \code{len_t} by \code{sz}\eqn{\times}{%\times}\code{sz}).
#' @details
#' Input \code{I} should has structure 'T_SS_mat', matrix with dimension
#' \code{len_t} by \code{sz}\eqn{\times}{%\times}\code{sz}.
#'
#' Input \code{pos} should be the position matrix with dimension
#' \code{M}\eqn{\times}{%\times}\code{len_t}. See \code{\link{bm_particle_intensity}},
#' \code{\link{ou_particle_intensity}}, \code{\link{fbm_particle_intensity}},
#' \code{\link{fbm_ou_particle_intensity}}.
#'
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
#' @keywords internal
fill_intensity <- function(len_t, M, I, pos, Ic, sz, sigma_p){
  for(i in 1:len_t){
    for(j in 1:M){
      xp = pos[j+M*(i-1),1]
      yp = pos[j+M*(i-1),2]
      x_range = floor(xp-3*sigma_p):ceiling(xp+3*sigma_p)
      y_range = floor(yp-3*sigma_p):ceiling(yp+3*sigma_p)
      x = rep(x_range,length(x_range))
      y = rep(y_range,each=length(x_range))
      dist_2 = (x-xp)^2+(y-yp)^2
      binary_result = (dist_2<=((3*sigma_p)^2))
      Ip = Ic[j]*exp(-dist_2 / (2*sigma_p^2))
      x_fill = x[binary_result]
      y_fill = y[binary_result]
      index_fill = y_fill+sz[1]*(x_fill-1)
      Ip_fill = Ip[binary_result]

      legitimate_index = (index_fill>0) & (index_fill<(sz[1]*sz[2]))
      if (length(legitimate_index) > 0){
        I[i,index_fill[legitimate_index]] = I[i,index_fill[legitimate_index]]+Ip_fill[legitimate_index]
      }
    }
  }
  return(I)
}

#' Compute numerical MSD
#' @description
#' Compute numerical mean squared displacement(MSD) based on particle trajectory.
#'
#'
#' @param pos position matrix for particle trajectory. See 'Details'.
#' @param M number of particles
#' @param len_t number of time steps
#'
#' @return A vector of numerical MSD for given lag times.
#' @details
#' Input \code{pos} should be the position matrix with dimension
#' \code{M}\eqn{\times}{%\times}\code{len_t}. See \code{\link{bm_particle_intensity}},
#' \code{\link{ou_particle_intensity}}, \code{\link{fbm_particle_intensity}},
#' \code{\link{fbm_ou_particle_intensity}}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' # Simulate particle trajectory for BM
#' M = 10
#' len_t = 50
#' sigma = 0.5
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = bm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma)
#'
#' # Compute numerical MSD
#' (num_msd = numerical_msd(pos=pos, M=M, len_t = len_t))
#' @keywords internal
numerical_msd <- function(pos, M,len_t){
  pos_msd = array(pos, dim=c(M, len_t, 2))
  msd_i = matrix(NaN,nrow=M,ncol=len_t-1)
  for(dt in 1:(len_t-1)){
    ndt = len_t-dt
    xdiff = pos_msd[,1:ndt,1]-pos_msd[,(1+dt):(ndt+dt),1]
    ydiff = pos_msd[,1:ndt,2]-pos_msd[,(1+dt):(ndt+dt),2]
    mean_square = xdiff^2+ydiff^2
    if (length(dim(mean_square))>1){
      msd_i[,dt] = apply(mean_square,1,function(x){mean(x,na.rm=T)})
    }else{msd_i[,dt] = mean_square}
  }
  #result_list = list()
  num_msd_mean = apply(msd_i,2,function(x){mean(x,na.rm=T)})
  #result_list$num_msd_mean = num_msd_mean
  #result_list$num_msd = msd_i
  num_msd_mean = c(0,num_msd_mean)
  return(num_msd_mean)
}

#' Compute anisotropic numerical MSD
#' @description
#' Compute numerical mean squared displacement(MSD) based on particle trajectory
#' for anisotropic processes in x,y-directions separately.
#'
#'
#' @param pos position matrix for particle trajectory. See 'Details'.
#' @param M number of particles
#' @param len_t number of time steps
#'
#' @return A matrix of numerical MSD for given lag times in x,y-directions,
#' dimension 2 by \code{len_t}.
#' @details
#' Input \code{pos} should be the position matrix with dimension
#' \code{M}\eqn{\times}{%\times}\code{len_t}. See \code{\link{bm_particle_intensity}},
#' \code{\link{ou_particle_intensity}}, \code{\link{fbm_particle_intensity}},
#' \code{\link{fbm_ou_particle_intensity}}.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' # Simulate particle trajectory for BM
#' M = 10
#' len_t = 50
#' sigma = c(0.5,0.1)
#' pos0 = matrix(100/8+0.75*100*runif(M*2),nrow=M,ncol=2)
#' pos = anisotropic_bm_particle_intensity(pos0=pos0,M=M,len_t=len_t,sigma=sigma)
#'
#' # Compute numerical MSD
#' (num_msd = anisotropic_numerical_msd(pos=pos, M=M, len_t=len_t))
#' @keywords internal
anisotropic_numerical_msd <- function(pos, M,len_t){
  pos_msd = array(pos, dim=c(M, len_t, 2))
  msd_x_i = matrix(NaN,nrow=M,ncol=len_t-1)
  msd_y_i = matrix(NaN,nrow=M,ncol=len_t-1)
  num_msd_mean = matrix(NaN, nrow=len_t,ncol=2)
  num_msd_mean[1,] = c(0,0)
  for(dt in 1:(len_t-1)){
    ndt = len_t-dt
    xdiff = pos_msd[,1:ndt,1]-pos_msd[,(1+dt):(ndt+dt),1]
    ydiff = pos_msd[,1:ndt,2]-pos_msd[,(1+dt):(ndt+dt),2]
    mean_square_x = xdiff^2
    mean_square_y = ydiff^2
    if (length(dim(mean_square_x))>1){
      msd_x_i[,dt] = apply(mean_square_x,1,function(x){mean(x,na.rm=T)})
    }else{msd_x_i[,dt] = mean_square_x}
    if (length(dim(mean_square_y))>1){
      msd_y_i[,dt] = apply(mean_square_y,1,function(x){mean(x,na.rm=T)})
    }else{msd_y_i[,dt] = mean_square_y}
    }
  #result_list = list()
  num_msd_mean[-1,1] = apply(msd_x_i,2,function(x){mean(x,na.rm=T)})
  num_msd_mean[-1,2] = apply(msd_y_i,2,function(x){mean(x,na.rm=T)})
  return(num_msd_mean)
}

#' Plot 2D particle trajectory
#' @description
#' Function to plot the particle trajectory after the \code{simulation} class
#' has been constructed.
#'
#' @param object an S4 object of class \code{simulation}
#' @param title main title of the plot. If \code{NA}, title is "model_name with
#' M particles" with \code{model_name} and \code{M} being field in \code{simulation}
#' class.
#'
#' @return 2D plot of particle trajectory for a given simulation from \code{simulation}
#' class.
#'
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#' plot_traj(sim_bm)
plot_traj<- function(object, title=NA){
  if(class(object)[1]=="simulation" || class(object)[1]=="aniso_simulation"){
    if(is.na(title)==T){
      if(class(object)[1]=="simulation"){title=paste(object@model_name,"with",object@M,"particles")
      }else if(class(object)[1]=="aniso_simulation"){
        title=paste("Anisotropic",object@model_name,"with",object@M,"particles")
      }
      }
    # highlight start and end point?
    traj1 = object@pos[seq(1,dim(object@pos)[1],by=object@M),]
    # change plot as 0,0 to be the top left corner
    #plot(traj1[,1],object@sz[1]-traj1[,2],ylim=c(0,object@sz[1]),
    plot(traj1[,1],traj1[,2],ylim=c(0,object@sz[1]),
         xlim=c(0,object@sz[2]),type="l",col=1,xlab="frame size",ylab="frame size",
         main=title)
    for (i in 2:object@M){
      v = object@pos[seq(i,dim(object@pos)[1],by=object@M),]
      lines(v[,1],v[,2],ylim=c(0,object@sz[1]),
            xlim=c(0,object@sz[2]),type="l", col=i)
    }
  }else{
    stop("Please input a simulation or aniso_simulation class object. \n")
  }
}

#' Show simulation object
#' @description
#' Function to print the \code{simulation} class object after the \code{simulation}
#' model has been constructed.
#'
#' @param object an S4 object of class \code{simulation}
#'
#' @return Show a list of important parameters in class \code{simulation}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#'
#' # Simulate simple diffusion for 100 images with 100 by 100 pixels
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
show.simulation <- function(object){
  cat("Frame size: ",object@sz, "\n")
  cat("Number of time steps: ",object@len_t, "\n")
  cat("Number of particles: ",object@M, "\n")
  cat("Stochastic process: ",object@model_name, "\n")
  cat("Variance of background noise: ",object@sigma_2_0, "\n")
  if(object@model_name == "BM"){
    cat("sigma_bm: ",object@param, "\n")
  }else if(object@model_name == "OU"){
    cat("(rho, sigma_ou): ",object@param,"\n")
  }else if(object@model_name == "FBM"){
    cat("(sigma_fbm, Hurst parameter): ",object@param, "\n")
  }else if(object@model_name == "OU+FBM"){
    cat("(rho, sigma_ou,sigma_fbm, Hurst parameter): ",object@param, "\n")
  }
}

#' Show anisotropic simulation object
#' @description
#' Function to print the \code{aniso_simulation} class object after the
#' \code{aniso_simulation} model has been constructed.
#'
#' @param object an S4 object of class \code{aniso_simulation}
#'
#' @return Show a list of important parameters in class \code{aniso_simulation}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#'
#' # Simulate simple diffusion for 100 images with 100 by 100 pixels
#' aniso_sim_bm = aniso_simulation(sz=100,len_t=100,sigma_bm=c(0.5,0.1))
#' show(aniso_sim_bm)
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
show.aniso_simulation <- function(object){
  cat("Frame size: ",object@sz, "\n")
  cat("Number of time steps: ",object@len_t, "\n")
  cat("Number of particles: ",object@M, "\n")
  cat("Stochastic process: ",object@model_name, "\n")
  cat("Variance of background noise: ",object@sigma_2_0, "\n")
  if(object@model_name == "BM"){
    cat("sigma_bm: ",object@param, "\n")
  }else if(object@model_name == "OU"){
    cat("(rho, sigma_ou): ",object@param,"\n")
  }else if(object@model_name == "FBM"){
    cat("(sigma_fbm, Hurst parameter): ",object@param, "\n")
  }else if(object@model_name == "OU+FBM"){
    cat("(rho, sigma_ou,sigma_fbm, Hurst parameter): ",object@param, "\n")
  }
}

#' Show scattering analysis of microscopy (SAM) object
#' @description
#' Function to print the \code{SAM} class object after the \code{SAM} model has
#' been constructed.
#'
#' @param object an S4 object of class \code{SAM}
#'
#' @return Show a list of important parameters in class \code{SAM}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#'
#' ## Simulate BM and get estimated parameters using BM model
#' # Simulation
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#'
#' # AIUQ method: fitting using BM model
#' sam = SAM(sim_object=sim_bm)
#' show(sam)
#'
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
show.sam <- function(object){
  cat("Fitted model: ",object@model_name, "\n")
  cat("Number of q ring: ",object@len_q, "\n")
  cat("Index of wave number selected: ",object@index_q, "\n")
  cat("True parameters in the model: ",object@param_truth, "\n")
  cat("Estimated parameters in the model: ",object@param_est, "\n")
  cat("True variance of background noise: ",object@sigma_2_0_truth, "\n")
  cat("Estimated variance of background noise: ",object@sigma_2_0_est, "\n")
  if(object@method=="AIUQ"){
    cat("Maximum log likelihood value: ",object@mle, "\n")
    cat("Akaike information criterion score: ",object@AIC, "\n")
  }
}

#' Show scattering analysis of microscopy for anisotropic processes (aniso_SAM) object
#' @description
#' Function to print the \code{aniso_SAM} class object after the
#' \code{aniso_SAM} model has been constructed.
#'
#' @param object an S4 object of class \code{aniso_SAM}
#'
#' @return Show a list of important parameters in class \code{aniso_SAM}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#'
#' ## Simulate BM and get estimated parameters using BM model
#' # Simulation
#' aniso_sim_bm = aniso_simulation(sz=100,len_t=100,sigma_bm=c(0.5,0.3))
#' show(aniso_sim_bm)
#'
#' # AIUQ method: fitting using BM model
#' aniso_sam = aniso_SAM(sim_object=aniso_sim_bm, AIUQ_thr=c(0.99,0))
#' show(aniso_sam)
#'
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
show.aniso_sam <- function(object){
  cat("Fitted model: ",object@model_name, "\n")
  cat("Number of q ring: ",object@len_q, "\n")
  cat("Index of wave number selected: ",object@index_q, "\n")
  cat("True parameters in the model: ",object@param_truth, "\n")
  cat("Estimated parameters in the model: ",object@param_est, "\n")
  cat("True variance of background noise: ",object@sigma_2_0_truth, "\n")
  cat("Estimated variance of background noise: ",object@sigma_2_0_est, "\n")
  if(object@method=="AIUQ"){
    cat("Maximum log likelihood value: ",object@mle, "\n")
    cat("Akaike information criterion score: ",object@AIC, "\n")
  }
}

#' Plot estimated MSD with uncertainty from SAM class
#' @description
#' Function to plot estimated MSD with uncertainty from \code{SAM} class, versus
#' true mean squared displacement(MSD) or given reference values.
#'
#' @param object an S4 object of class \code{SAM}
#' @param msd_truth  a vector/matrix of true MSD or reference MSD value,
#' default is \code{NA}
#' @param title main title of the plot. If \code{NA}, title is "model_name" with
#' \code{model_name} being a field in \code{SAM} class representing fitted model.
#' @param log10 a logical evaluating to TRUE or FALSE indicating whether a plot
#' in log10 scale is generated
#'
#' @return A plot of estimated MSD with uncertainty versus truth/reference values.
#' @export
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#'
#' ## Simulate BM and get estimated parameters with uncertainty using BM model
#' # Simulation
#' set.seed(1)
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#'
#' # AIUQ method: fitting using BM model
#' sam = SAM(sim_object=sim_bm, uncertainty=TRUE,AIUQ_thr=c(0.999,0))
#' show(sam)
#'
#' plot_MSD(object=sam, msd_truth=sam@msd_truth) #in log10 scale
#' plot_MSD(object=sam, msd_truth=sam@msd_truth,log10=FALSE) #in real scale
#'
#' @references
#' Gu, M., He, Y., Liu, X., & Luo, Y. (2023). Ab initio uncertainty
#' quantification in scattering analysis of microscopy.
#' arXiv preprint arXiv:2309.02468.
#'
#' Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021).
#' Uncertainty quantification and estimation in differential dynamic microscopy.
#' Physical Review E, 104(3), 034610.
plot_MSD<-function(object, msd_truth=NA, title=NA, log10=TRUE){
  if(class(object)[1]=='SAM'){
    if(is.na(title)==T){title=paste(object@model_name)}
    if(!is.na(msd_truth)[1]){
      len_t = min(length(msd_truth),length(object@d_input))
      msd_truth_here = msd_truth[2:len_t]
      if(log10==T){
        plot(log10(object@d_input[2:len_t]),log10(msd_truth_here),
             type='l',col='black',ylab='log10(MSD)',
             xlab=expression(paste("log10(", Delta, "t)",sep="")), main=title,
             lwd=2,lty=1, ylim = c(log10(min(msd_truth_here)),log10(max(msd_truth_here))*1.1),
             mgp=c(2.5,1,0))
        if(!is.na(object@msd_upper)[1]){
          polygon(log10(c(object@d_input[-1],rev(object@d_input[-1]))),
                  log10(c(object@msd_upper[-1],rev(object@msd_lower[-1]))),
                  col = "grey80", border = F)
        }
        lines(log10(object@d_input[-1]),log10(object@msd_est[-1]),type='p',
              col='blue', pch=20, cex=0.8)
        lines(log10(object@d_input[2:len_t]),log10(msd_truth_here),
              type='l',col='black',lwd=2)
        legend('topleft',lty=c(1,NA),pch=c(NA,20),col=c('black','blue'),
               legend=c('Reference','SAM'),lwd=c(2,2), cex=0.7)
      }else{
        plot(object@d_input[2:len_t],msd_truth_here,type='l',col='black',
             ylab='MSD',xlab=expression(Delta~t),main=title,lwd=2,lty=1,
             ylim = c(0,max(msd_truth_here)*1.2),mgp=c(2.5,1,0))
        if(!is.na(object@msd_upper)[1]){
          polygon(c(object@d_input[-1],rev(object@d_input[-1])),
                  c(object@msd_upper[-1],rev(object@msd_lower[-1])),
                  col = "grey80", border = F)
        }
        lines(object@d_input[-1],object@msd_est[-1],type='p',col='blue', pch=20, cex=0.8)
        lines(object@d_input[2:len_t],msd_truth_here,type='l',col='black',lwd=2)
        legend('topleft',lty=c(1,NA),pch=c(NA,20), col=c('black','blue'),
               legend=c('Reference','SAM'),lwd=c(2,2), cex=0.7)
      }
    }else{
      if(log10==T){
        plot(log10(object@d_input[-1]),log10(object@msd_est[-1]),type='p',
             col='blue', pch=20, cex=0.8,ylab='log10(MSD)',
             xlab=expression(paste("log10(", Delta, "t)",sep="")),main=title,
             lwd=2,lty=1, ylim = c(log10(min(object@msd_est[-1])),log10(max(object@msd_est))*1.2),
             mgp=c(2.5,1,0))
        if(!is.na(object@msd_upper)[1]){
          polygon(log10(c(object@d_input[-1],rev(object@d_input[-1]))),
                  log10(c(object@msd_upper[-1],rev(object@msd_lower[-1]))),
                  col = "grey80", border = F)
        }
        lines(log10(object@d_input[-1]),log10(object@msd_est[-1]),type='p',col='blue', pch=20, cex=0.8)
        legend('topleft',pch=20, lty=NA,col=c('blue'),
               legend=c('SAM'),lwd=c(2), cex=0.7)
      }else{
        plot(object@d_input[-1],object@msd_est[-1],type='p',col='blue', pch=20, cex=0.8,
             ylab='MSD',xlab=expression(Delta~t),main=title,lwd=2,lty=1,
             ylim = c(0,max(object@msd_est[-1])*1.1),mgp=c(2.5,1,0))
        if(!is.na(object@msd_upper)[1]){
          polygon(c(object@d_input[-1],rev(object@d_input[-1])),
                  c(object@msd_upper[-1],rev(object@msd_lower[-1])),col = "grey80", border = F)
        }
        lines(object@d_input[-1],object@msd_est[-1],type='p',col='blue', pch=20, cex=0.8)
        legend('topleft',pch=20,lty=NA, col=c('blue'),legend=c('SAM'),
               lwd=c(2), cex=0.7)
      }
    }
  }else if(class(object)[1]=='aniso_SAM'){
    if(is.na(title)==T){title=paste("Anisotropic",object@model_name)}
    if(!is.na(msd_truth)[1]){
      if(ncol(msd_truth)!=2){stop("Please input a matrix with 2 columns where
                                  each column holds MSD truth for x,y directions,
                                  respectively.\n")}
      len_t = min(nrow(msd_truth),length(object@d_input))
      msd_true_x = msd_truth[2:len_t,1]
      msd_true_y = msd_truth[2:len_t,2]
      msd_x = object@msd_est[-1,1]
      msd_y = object@msd_est[-1,2]
      if(log10==T){
        plot(log10(object@d_input[2:len_t]),log10(msd_true_x),
             type='l',col='black',ylab='log10(MSD)',
             xlab=expression(paste("log10(", Delta, "t)",sep="")), main=title,
             lwd=2,lty=1, ylim = c(log10(min(msd_true_x,msd_true_y)),
                                   log10(max(msd_true_x,msd_true_y))*1.1),
             mgp=c(2.5,1,0))
        lines(log10(object@d_input[2:len_t]),log10(msd_true_y),col='black',lwd=2,lty=2)
        if(!is.na(object@msd_x_upper)[1]){
          polygon(log10(c(object@d_input[-1],rev(object@d_input[-1]))),
                  log10(c(object@msd_x_upper[-1],rev(object@msd_x_lower[-1]))),
                  col = "grey80", border = F)
        }
        if(!is.na(object@msd_y_upper)[1]){
          polygon(log10(c(object@d_input[-1],rev(object@d_input[-1]))),
                  log10(c(object@msd_y_upper[-1],rev(object@msd_y_lower[-1]))),
                  col = "grey80", border = F)
        }
        lines(log10(object@d_input[-1]),log10(msd_x),type='p',col='blue', pch=20, cex=0.8)
        lines(log10(object@d_input[-1]),log10(msd_y),type='p',col='blue', pch=17, cex=0.8)
        lines(log10(object@d_input[2:len_t]),log10(msd_true_x),col='black',lwd=2,lty=1)
        lines(log10(object@d_input[2:len_t]),log10(msd_true_y),col='black',lwd=2,lty=2)
        legend('topleft',lty=c(1,2,NA,NA),pch=c(NA,NA,20,17),
               col=c('black','black','blue','blue'),
               legend=c('Reference x','Reference y','SAM x','SAM y'),lwd=c(2,2,2,2),
               cex=0.7)

      }else{
        plot(object@d_input[2:len_t],msd_true_x,
             type='l',col='black',ylab='MSD',
             xlab=expression(paste("", Delta, "t",sep="")), main=title,
             lwd=2,lty=1, ylim = c(min(msd_true_x,msd_true_y),
                                   max(msd_true_x,msd_true_y)*1.1),
             mgp=c(2.5,1,0))
        lines(object@d_input[2:len_t],msd_true_y,col='black',lwd=2,lty=2)
        if(!is.na(object@msd_x_upper)[1]){
          polygon(c(object@d_input[-1],rev(object@d_input[-1])),
                  c(object@msd_x_upper[-1],rev(object@msd_x_lower[-1])),
                  col = "grey80", border = F)
        }
        if(!is.na(object@msd_y_upper)[1]){
          polygon(c(object@d_input[-1],rev(object@d_input[-1])),
                  c(object@msd_y_upper[-1],rev(object@msd_y_lower[-1])),
                  col = "grey80", border = F)
        }
        lines(object@d_input[-1],msd_x,type='p',col='blue', pch=20, cex=0.8)
        lines(object@d_input[-1],msd_y,type='p',col='blue', pch=17, cex=0.8)
        lines(object@d_input[2:len_t],msd_true_x,col='black',lwd=2,lty=1)
        lines(object@d_input[2:len_t],msd_true_y,col='black',lwd=2,lty=2)
        legend('topleft',lty=c(1,2,NA,NA),pch=c(NA,NA,20,17),
               col=c('black','black','blue','blue'),
               legend=c('Reference x','Reference y','SAM x','SAM y'),lwd=c(2,2,2,2),
               cex=0.7)

      }
    }else{
      if(log10==T){
        msd_x = object@msd_est[-1,1]
        msd_y = object@msd_est[-1,2]
        plot(log10(object@d_input[-1]),log10(msd_x),type='p', col='blue',pch=20,
             cex=0.8, main=title, lwd=2, lty=1, ylab='log10(MSD)',
             xlab=expression(paste("log10(", Delta, "t)",sep="")),
             ylim = c(log10(min(msd_x,msd_y)),log10(max(msd_x,msd_y))*1.2),
             mgp=c(2.5,1,0))
        lines(log10(object@d_input[-1]),log10(msd_y),type='p',col='blue', pch=17, cex=0.8)
        if(!is.na(object@msd_x_upper)[1]){
          polygon(log10(c(object@d_input[-1],rev(object@d_input[-1]))),
                  log10(c(object@msd_x_upper[-1],rev(object@msd_x_lower[-1]))),
                  col = "grey80", border = F)}
        if(!is.na(object@msd_y_upper)[1]){
          polygon(log10(c(object@d_input[-1],rev(object@d_input[-1]))),
                  log10(c(object@msd_y_upper[-1],rev(object@msd_y_lower[-1]))),
                  col = "grey80", border = F)}
        lines(log10(object@d_input[-1]),log10(msd_x),type='p',col='blue', pch=20, cex=0.8)
        lines(log10(object@d_input[-1]),log10(msd_y),type='p',col='blue', pch=17, cex=0.8)

        legend('topleft',pch=c(20,17), lty=c(NA,NA),col=c('blue','blue'),
               legend=c('SAM x','SAM y'),lwd=c(2,2), cex=0.7)
      }else{
        msd_x = object@msd_est[-1,1]
        msd_y = object@msd_est[-1,2]
        plot(object@d_input[-1],msd_x,type='p', col='blue',pch=20,
             cex=0.8, main=title, lwd=2, lty=1, ylab='MSD',
             xlab=expression(paste( Delta, "t",sep="")),
             ylim = c(min(msd_x,msd_y),max(msd_x,msd_y)*1.2),
             mgp=c(2.5,1,0))
        lines(object@d_input[-1],msd_y,type='p',col='blue', pch=17, cex=0.8)
        if(!is.na(object@msd_x_upper)[1]){
          polygon(c(object@d_input[-1],rev(object@d_input[-1])),
                  c(object@msd_x_upper[-1],rev(object@msd_x_lower[-1])),
                  col = "grey80", border = F)}
        if(!is.na(object@msd_y_upper)[1]){
          polygon(c(object@d_input[-1],rev(object@d_input[-1])),
                  c(object@msd_y_upper[-1],rev(object@msd_y_lower[-1])),
                  col = "grey80", border = F)}
        lines(object@d_input[-1],msd_x,type='p',col='blue', pch=20, cex=0.8)
        lines(object@d_input[-1],msd_y,type='p',col='blue', pch=17, cex=0.8)

        legend('topleft',pch=c(20,17), lty=c(NA,NA),col=c('blue','blue'),
               legend=c('SAM x','SAM y'),lwd=c(2,2), cex=0.7)

      }
    }
  }
  else{stop("Please input an SAM or aniso_SAM class object. \n")}
}


#' Plot 2D intensity
#' @description
#' Function to plot 2D intensity profile for a certain frame, default is to plot
#' the first frame. Input can be a matrix (2D) or an array (3D).
#'
#' @param intensity intensity profile
#' @param intensity_str structure of the intensity profile, options from
#' ('SST_array','S_ST_mat','T_SS_mat', 'SS_T_mat'). See 'Details'.
#' @param frame frame index
#' @param title main title of the plot. If \code{NA}, title is "intensity profile
#' for frame n" with n being the frame index in \code{frame}.
#' @param color a logical evaluating to TRUE or FALSE indicating whether a colorful
#' plot is generated
#' @param sz frame size of simulated image with default \code{c(200,200)}.
#'
#' @return 2D plot in gray scale (or with color) of selected frame.
#' @details
#' By default \code{intensity_str} is set to 'T_SS_mat', a time by space\eqn{\times}{%\times}space
#' matrix, which is the structure of intensity profile obtained from \code{simulation}
#' class. For \code{intensity_str='SST_array'} , input intensity profile should be a
#' space by space by time array, which is the structure from loading a tif file.
#' For \code{intensity_str='S_ST_mat'}, input intensity profile should be a
#' space by space\eqn{\times}{%\times}time matrix. For \code{intensity_str='SS_T_mat'},
#' input intensity profile should be a space\eqn{\times}{%\times}space by time matrix.
#'
#' @author \packageAuthor{AIUQ}
#' @examples
#' library(AIUQ)
#' sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
#' show(sim_bm)
#' plot_intensity(sim_bm@intensity, sz=sim_bm@sz)
#'
#' @export
plot_intensity<-function(intensity,intensity_str="T_SS_mat",frame=1,sz=NA,
                         title=NA, color=FALSE){
  if(is.na(title)){
    title = paste("intensity profile for frame ", frame, sep="")
  }
  if(length(sz)==1 && is.na(sz)){
    intensity_list = intensity_format_transform(intensity=intensity,
                                                intensity_str=intensity_str)
    intensity_trans = intensity_list$intensity
    sz_x = intensity_list$sz_x
    sz_y = intensity_list$sz_y
  }else{
    intensity_list = intensity_format_transform(intensity=intensity,
                                                intensity_str=intensity_str,sz=sz)
    intensity_trans = intensity_list$intensity
    sz_x = intensity_list$sz_x
    sz_y = intensity_list$sz_y
  }
  if(color==FALSE){
    plot_m = t(matrix(intensity_trans[,frame],sz_y,sz_x))
    #reversed_m =  plot_m[, ncol(plot_m):1]
    plot3D::image2D(plot_m,main=title,col = grey(seq(0, 1, length = 256)))
  }else{
    plot_m = t(matrix(intensity_trans[,frame],sz_y,sz_x))
    #reversed_m =  plot_m[, ncol(plot_m):1]
    plot3D::image2D(plot_m,main=title)
  }
}


#' Compute dynamic image structure function
#' @description
#' Compute dynamic image structure function(Dqt) using Fourier transformed
#' intensity profile and a selection of wave number(q) range.
#'
#' @param len_q number of wave number
#' @param index_q a vector of selected wave number index
#' @param len_t number of time steps
#' @param I_q_matrix intensity profile in reciprocal space (after Fourier transformation)
#' @param q_ori_ring_loc_unique_index index for wave vector that give unique frequency
#' @param sz frame size of intensity profile
#'
#' @return Matrix of dynamic image structure with dimension \code{len_q} by \code{len_t-1}.
#' @export
#' @author \packageAuthor{AIUQ}
#' @details
#' Dynamic image structure function(Dqt) can be obtained from ensemble average
#' of absolute values squared of Four transformed intensity difference:
#' \deqn{D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}{%D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}
#' See 'References'.
#'
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
#' @keywords internal
SAM_Dqt<-function(len_q,index_q,len_t,I_q_matrix,q_ori_ring_loc_unique_index,sz){
  Dqt = matrix(NA,len_q,len_t-1)
  for (q_j in index_q){
    I_q_cur = I_q_matrix[q_ori_ring_loc_unique_index[[q_j]],]
    for (t_i in 1:(len_t-1)){
      Dqt[q_j,t_i]=mean((abs(I_q_cur[,(t_i+1):len_t]-I_q_cur[,1:(len_t-t_i)]))^2/(sz[1]*sz[2]),na.rm=T)
    }
  }
  return(Dqt)
}

#' Compute l2 loss for Dqt with fixed A(q) and B
#' @description
#' Compute l2 loss for dynamic image structure function(Dqt) using fixed A(q)
#' and B parameters.
#'
#' @param param a vector of natural logarithm of parameters
#' @param Dqt_cur observed dynamic image structure function. See 'Details'.
#' @param q_cur wave vector in unit of um^-1
#' @param A_est_q_cur estimated value of A(q). This parameter is determined by
#' the properties of the imaged material and imaging optics. See 'References'.
#' @param B_est estimated value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#' @param msd_fn msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return Squared differences between the true Dqt and the predicted Dqt.
#' @details
#' Dynamic image structure function(Dqt) can be obtained from ensemble average
#' of absolute values squared of Four transformed intensity difference:
#' \deqn{D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}{%D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}
#' See 'References'.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @keywords internal
l2_fixedAB<-function(param,Dqt_cur,q_cur,A_est_q_cur,B_est,
                     d_input,model_name,msd_fn=NA,msd_grad_fn=NA){
  theta = exp(param)
  msd_list = get_MSD_with_grad(theta=theta,d_input=d_input[-1],
                               model_name=model_name,msd_fn,msd_grad_fn=msd_grad_fn)
  msd = msd_list$msd
  sum((Dqt_cur-(A_est_q_cur*(1-exp(-q_cur^2*msd/4))+B_est))^2)
}

#' Compute l2 loss for Dqt
#' @description
#' Compute l2 loss for dynamic image structure function(Dqt) using A(q) and B
#' are both estimated within the model.
#'
#' @param param a vector of natural logarithm of parameters
#' @param Dqt_cur observed dynamic image structure function. See 'Details'.
#' @param q_cur wave vector in unit of um^-1
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#' @param msd_fn msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return Squared differences between the true Dqt and the predicted Dqt.
#' @details
#' Dynamic image structure function(Dqt) can be obtained from ensemble average
#' of absolute values squared of Four transformed intensity difference:
#' \deqn{D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}{%D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}
#' See 'References'.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @keywords internal
l2_estAB<-function(param,Dqt_cur,q_cur,d_input,model_name,
                   msd_fn=NA,msd_grad_fn=NA){
  theta = exp(param)
  A_cur = theta[length(theta)]
  B_cur= theta[length(theta)-1]
  theta_msd = theta[-((length(theta)-1):length(theta))]
  msd_list = get_MSD_with_grad(theta=theta_msd,d_input=d_input[-1],
                               model_name=model_name,msd_fn,msd_grad_fn=msd_grad_fn)
  msd = msd_list$msd
  sum((Dqt_cur-(A_cur*(1-exp(-q_cur^2*msd/4))+B_cur))^2)
}

#' Minimize l2 loss for Dqt with fixed A(q) and B
#' @description
#' Minimize l2 loss function for dynamic image structure function(Dqt) with
#' fixed A(q) and B, and return estimated parameters and mean squared
#' displacement(MSD).
#'
#' @param param a vector of natural logarithm of parameters
#' @param q wave vector in unit of um^-1
#' @param index_q selected index of wave number
#' @param Dqt observed dynamic image structure function. See 'Details'.
#' @param A_est_q estimated value of A(q). This parameter is determined by
#' the properties of the imaged material and imaging optics. See 'References'.
#' @param B_est estimated value of B. This parameter is determined by the noise
#' in the system. See 'References'.
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#' @param msd_fn msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return A list of estimated parameters and MSD from minimizing the l2 loss
#' function.
#' @details
#' Dynamic image structure function(Dqt) can be obtained from ensemble average
#' of absolute values squared of Four transformed intensity difference:
#' \deqn{D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}{%D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}
#' See 'References'.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @keywords internal
theta_est_l2_dqt_fixedAB<-function(param,q,index_q,Dqt,A_est_q,B_est,
                                   d_input,model_name,msd_fn=NA,msd_grad_fn=NA){
  param_est_l2 = matrix(nrow=length(q),ncol=length(param))
  for(q_j in index_q){
    m_optim=try(optim(param,l2_fixedAB,Dqt_cur=Dqt[q_j,],d_input=d_input,
                      A_est_q_cur=A_est_q[q_j],B_est=B_est,q_cur=q[q_j],
                      model_name=model_name,msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,
                      method='L-BFGS-B' ),silent=T)
    try_num=0
    while(!is.numeric(m_optim[[1]])){
      try_num=try_num+1
      param=param+runif(length(param))
      m_optim=try(optim(param,l2_fixedAB,Dqt_cur=Dqt[q_j,],d_input=d_input,
                        A_est_q_cur=A_est_q[q_j],B_est=B_est,q_cur=q[q_j],
                        model_name=model_name,msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,
                        method='L-BFGS-B' ),silent=T)
      if(try_num>10){
        break
      }
    }
    ##try no gradient if still no value
    try_num=0
    while(!is.numeric(m_optim[[1]])){
      try_num=try_num+1
      param=param+runif(length(param))
      m_optim=try(optim(param,l2_fixedAB,Dqt_cur=Dqt[q_j,],d_input=d_input,
                        A_est_q_cur=A_est_q[q_j],B_est=B_est,q_cur=q[q_j],
                        model_name=model_name,msd_fn=msd_fn,msd_grad_fn=msd_grad_fn),silent=T)
      if(try_num>10){
        break
      }
    }
    param_est_l2[q_j,]=(m_optim$par)
  }
  param_ddm = apply(exp(param_est_l2),1,function(x){get_est_param(x,model_name=model_name)})
  if(is.vector(param_ddm)==T){
    param_ddm = mean(param_ddm,na.rm=T)
  }else{
    param_ddm = apply(param_ddm,1,function(x){mean(x,na.rm=T)})
  }
  msd_ddm = get_MSD(theta=param_ddm,d_input=d_input,model_name=model_name,msd_fn=msd_fn)
  ddm_result = list()
  ddm_result$param_est = param_ddm
  ddm_result$msd_est = msd_ddm
  return(ddm_result)
}

#' Minimize l2 loss for Dqt
#' @description
#' Minimize l2 loss function for dynamic image structure function(Dqt), and
#' return estimated parameters and mean squared displacement(MSD).
#'
#' @param param a vector of natural logarithm of parameters
#' @param A_ini initial value of A(q) to be optimized over. Note true A(q) is
#' determined by the properties of the imaged material and imaging optics.
#' See 'References'.
#' @param q wave vector in unit of um^-1
#' @param index_q selected index of wave number
#' @param Dqt observed dynamic image structure function. See 'Details'.
#' @param d_input sequence of lag times
#' @param model_name model name for the fitted model, options from ('BM','OU',
#' 'FBM',OU+FBM','user_defined')
#' @param msd_fn msd_fn user defined mean squared displacement structure (MSD), a
#' function of \code{param} parameters and \code{d_input} lag times
#' @param msd_grad_fn user defined MSD gradient structure,  a function of
#' \code{param} and \code{d_input}
#'
#' @return A list of estimated parameters and MSD from minimizing the l2 loss
#' function.
#' @details
#' Dynamic image structure function(Dqt) can be obtained from ensemble average
#' of absolute values squared of Four transformed intensity difference:
#' \deqn{D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}{%D(q,\Delta t) = \langle |\Delta \hat{I}(q,t,\Delta t)|^2\rangle}
#' See 'References'.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @keywords internal
theta_est_l2_dqt_estAB<- function(param,A_ini,q,index_q,Dqt,d_input,
                                  model_name,msd_fn=NA,msd_grad_fn=NA){

  param_est_l2 = matrix(nrow=length(q),ncol=length(param)+1)
  for(q_j in index_q){
    param_start=c(param,log(abs(A_ini[q_j])))
    if(q_j==length(q)){
      param_start=c(param,0)
    }
    m_optim=try(optim(param_start,l2_estAB,Dqt_cur=Dqt[q_j,],d_input=d_input,
                      q_cur=q[q_j],model_name=model_name,msd_fn=msd_fn,
                      msd_grad_fn=msd_grad_fn,method='L-BFGS-B' ),silent=T)

    try_num=0
    while(!is.numeric(m_optim[[1]])){
      try_num=try_num+1
      param_start=param_start+runif(length(param)+1)
      m_optim=try(optim(param_start,l2_estAB,Dqt_cur=Dqt[q_j,],d_input=d_input,
                        q_cur=q[q_j],model_name=model_name,msd_fn=msd_fn,
                        msd_grad_fn=msd_grad_fn,method='L-BFGS-B' ),silent=T)
      if(try_num>10){
        break
      }
    }

    try_num=0
    while(!is.numeric(m_optim[[1]])){
      try_num=try_num+1
      param_start=param_start+runif(length(param)+1)
      m_optim=try(optim(param_start,l2_estAB,Dqt_cur=Dqt[q_j,],d_input=d_input,
                        q_cur=q[q_j],model_name=model_name,msd_fn=msd_fn,
                        msd_grad_fn=msd_grad_fn),silent=T)

      if(try_num>10){
        break
      }
    }

    param_est_l2[q_j,]=(m_optim$par)
  }

  param_ddm = apply(exp(param_est_l2[,-((length(param)):length(param)+1)]),1,
                    function(x){get_est_param(x,model_name=model_name)})
  if(is.vector(param_ddm)==T){
    param_ddm = mean(param_ddm,na.rm=T)
  }else{
    param_ddm = apply(param_ddm,1,function(x){mean(x,na.rm=T)})
  }
  msd_ddm = get_MSD(theta=param_ddm,d_input=d_input,model_name=model_name,msd_fn=msd_fn)
  ddm_result = list()
  ddm_result$param_est = param_ddm
  ddm_result$msd_est = msd_ddm
  ddm_result$sigma_2_0_est = mean(exp(param_est_l2[,length(param)]),na.rm=T)
  ddm_result$A_est=exp(param_est_l2[,length(param)+1])

  return(ddm_result)
}

#' Compute observed dynamic image structure function
#' @description
#' Compute observed dynamic image structure function (Dqt) using object of
#' \code{SAM} class.
#'
#' @param object an S4 object of class \code{SAM}
#' @param index_q wavevector range used for computing Dqt
#'
#' @return  A matrix of observed dynamic image structure function with dimension
#' \code{len_q} by \code{len_t-1}.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @examples
#' library(AIUQ)
#' sim_bm = simulation(len_t=100,sz=100,sigma_bm=0.5)
#' show(sim_bm)
#' sam = SAM(sim_object = sim_bm)
#' show(sam)
#' Dqt = get_dqt(object=sam)
get_dqt <- function(object, index_q = NA){
  if(class(object)[1]=='SAM'){
    len_q = object@len_q
    len_t = object@len_t
    sz = object@sz
    if(length(index_q)==1 && is.na(index_q)){
      index_q = 1:len_q
    }
    Dqt = matrix(NA,len_q,len_t-1)
    for (q_j in index_q){
      index_cur = object@q_ori_ring_loc_unique_index[[q_j]]
      I_q_cur = object@I_q[index_cur,]
      for (t_i in 1:(len_t-1)){
        Dqt[q_j,t_i]=mean((abs(I_q_cur[,(t_i+1):len_t]-I_q_cur[,1:(len_t-t_i)]))^2/(sz[1]*sz[2]),na.rm=T)
      }
    }
    return(Dqt)
  }else{stop("Please input an SAM class. \n")}
}

#' Compute empirical intermediate scattering function
#' @description
#' Compute empirical intermediate scattering function (ISF) using object of
#' \code{SAM} class.
#'
#' @param object an S4 object of class \code{SAM}
#' @param index_q wavevector range used for computing ISF
#'
#' @return  A matrix of empirical intermediate scattering function with dimension
#' \code{len_q} by \code{len_t-1}.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @examples
#' library(AIUQ)
#' sim_bm = simulation(len_t=100,sz=100,sigma_bm=0.5)
#' show(sim_bm)
#' sam = SAM(sim_object = sim_bm)
#' show(sam)
#' ISF = get_isf(object=sam)
get_isf <- function(object, index_q = NA){
  if(class(object)[1]=='SAM'){
    len_q = object@len_q
    len_t = object@len_t
    A_est = object@A_est_ini
    B_est = object@B_est_ini
    if(length(index_q)==1 && is.na(index_q)){
      index_q = 1:len_q
    }

    if(nrow(object@Dqt)*ncol(object@Dqt)==1 && is.na(object@Dqt)){
      sz = object@sz
      Dqt = matrix(NA,len_q,len_t-1)
      isf = matrix(NA,len_q,len_t-1)
      for (q_j in index_q){
        index_cur = object@q_ori_ring_loc_unique_index[[q_j]]
        I_q_cur = object@I_q[index_cur,]
        for (t_i in 1:(len_t-1)){
          Dqt[q_j,t_i]=mean((abs(I_q_cur[,(t_i+1):len_t]-I_q_cur[,1:(len_t-t_i)]))^2/(sz[1]*sz[2]),na.rm=T)
        }
        if(A_est[q_j]==0){break}
        isf[q_j,] = 1-(Dqt[q_j,]-B_est)/A_est[q_j]
      }
    }else{
      Dqt = object@Dqt
      isf = matrix(NA,len_q,len_t-1)
      for (q_j in index_q){
        if(A_est[q_j]==0){break}
        isf[q_j,] = 1-(Dqt[q_j,]-B_est)/A_est[q_j]
      }
    }
    return(isf)
  }else{stop("Please input an SAM class. \n")}
}

#' Compute modeled intermediate scattering function
#' @description
#' Compute modeled intermediate scattering function (ISF) using object of
#' \code{SAM} class.
#'
#' @param object an S4 object of class \code{SAM}
#' @param index_q wavevector range used for computing ISF
#'
#' @return  A matrix of modeled intermediate scattering function with dimension
#' \code{len_q} by \code{len_t-1}.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @examples
#' library(AIUQ)
#' sim_bm = simulation(len_t=100,sz=100,sigma_bm=0.5)
#' show(sim_bm)
#' sam = SAM(sim_object = sim_bm)
#' show(sam)
#' modeled_ISF = modeled_isf(object=sam)
modeled_isf <- function(object, index_q=NA){
  if(class(object)[1]=='SAM'){
    len_q = object@len_q
    len_t = object@len_t
    q = object@q
    msd_est = object@msd_est[-1]
    isf = matrix(NA,len_q,len_t-1)
    if(length(index_q)==1 && is.na(index_q)){
      index_q = 1:len_q
    }
    for(q_j in index_q){
      q_selected = q[q_j]
      isf[q_j,] = exp(-q_selected^2*msd_est/4)
    }
    return(isf)
  }else{stop("Please input an SAM class. \n")}
}

#' Compute modeled dynamic image structure function
#' @description
#' Compute modeled dynamic image structure function (Dqt) using object of
#' \code{SAM} class.
#'
#' @param object an S4 object of class \code{SAM}
#' @param index_q wavevector range used for computing Dqt
#'
#' @return  A matrix of modeled dynamic image structure function with dimension
#' \code{len_q} by \code{len_t-1}.
#'
#' @author \packageAuthor{AIUQ}
#' @export
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
#' @examples
#' library(AIUQ)
#' sim_bm = simulation(len_t=100,sz=100,sigma_bm=0.5)
#' show(sim_bm)
#' sam = SAM(sim_object = sim_bm)
#' show(sam)
#' modeled_Dqt = modeled_dqt(object=sam)
modeled_dqt <- function(object, index_q = NA){
  if(class(object)[1]=='SAM'){
    len_q = object@len_q
    len_t = object@len_t
    A_est = object@A_est_ini
    B_est = object@sigma_2_0_est*2
    q = object@q
    msd_est = object@msd_est[-1]
    if(length(index_q)==1 && is.na(index_q)){
      index_q = 1:len_q
    }
    if(nrow(object@modeled_ISF)*ncol(object@modeled_ISF)==1 && is.na(object@modeled_ISF)){
      isf = matrix(NA,len_q,len_t-1)
      dqt = matrix(NA,len_q,len_t-1)
      for(q_j in index_q){
        q_selected = q[q_j]
        isf[q_j,] = exp(-q_selected^2*msd_est/4)
        if(A_est[q_j]==0){break}
        dqt[q_j,] = A_est[q_j]*(1-isf[q_j,])+B_est
      }
    }else{
      isf = object@modeled_ISF
      dqt = matrix(NA,len_q,len_t-1)
      for (q_j in index_q){
        if(A_est[q_j]==0){break}
        dqt[q_j,] = A_est[q_j]*(1-isf[q_j,])+B_est
      }
    }
    return(dqt)
  }else{stop("Please input an SAM class. \n")}
}

