#' Scattering analysis of microscopy for anisotropic processes
#'
#' @description
#' Fast parameter estimation in scattering analysis of microscopy for anisotropic
#' processes, using AIUQ method.
#'
#' @param intensity intensity profile. See 'Details'.
#' @param intensity_str structure of the intensity profile, options from
#' ('SST_array','S_ST_mat','T_SS_mat'). See 'Details'.
#' @param sz frame size of the intensity profile in x and y directions,
#' number of pixels contained in each frame equals sz_x by sz_y.
#' @param pxsz size of one pixel in unit of micron, 1 for simulated data
#' @param mindt minimum lag time, 1 for simulated data
#' @param AIUQ_thr threshold for wave number selection, numeric vector of two
#' elements with values between 0 and 1. See 'Details'.
#' @param model_name fitted model, options from  ('BM','OU','FBM','OU+FBM',
#' 'user_defined'), with Brownian motion as the default model. See 'Details'.
#' @param sigma_0_2_ini initial value for background noise. If NA, use minimum
#' value of absolute square of intensity profile in reciprocal space.
#' @param msd_fn user defined mean squared displacement(MSD) structure, a
#' function of parameters and lag times. NA if \code{model_name} is not
#' 'user_defined'.
#' @param msd_grad_fn gradient for user defined mean squared displacement
#' structure. If \code{NA}, then numerical gradient will be used for parameter
#' estimation in \code{'user_defined'} model.
#' @param num_param number of parameters need to be estimated in the intermediate
#' scattering function, need to be non-NA value for user_defined' model.
#' @param param_initial initial values for param estimation.
#' @param num_optim number of optimization.
#' @param uncertainty a logical evaluating to TRUE or FALSE indicating whether
#' parameter uncertainty should be computed.
#' @param M number of particles. See 'Details'.
#' @param sim_object NA or an S4 object of class \code{simulation}.
#' @param msd_truth true MSD or reference MSD value.
#' @param method methods for parameter estimation, options from ('AIUQ', 'DDM').
#' @param index_q_AIUQ index range for wave number when using AIUQ method. See 'Details'.
#' @param message_out a logical evaluating to TRUE or FALSE indicating whether
#' or not to output the message.
#' @param square a logical evaluating to TRUE or FALSE indicating whether or not
#' to crop the original intensity profile into square image.
#'
#' @details
#' For simulated data using \code{aniso_simulation} in AIUQ package, \code{intensity}
#' will be automatically extracted from \code{aniso_simulation} class.
#'
#' By default \code{intensity_str} is set to 'T_SS_mat', a time by space\eqn{\times}{%\times}space
#' matrix, which is the structure of intensity profile obtained from \code{aniso_simulation}
#' class. For \code{intensity_str='SST_array'} , input intensity profile should be a
#' space by space by time array, which is the structure from loading a tif file.
#' For \code{intensity_str='S_ST_mat'}, input intensity profile should be a
#' space by space\eqn{\times}{%\times}time matrix.
#'
#' By default \code{AIUQ_thr} is set to \code{c(1,1)}, uses information from all
#' complete q rings. The first element affects maximum wave number selected,
#' and second element controls minimum proportion of wave number selected. By
#' setting 1 for the second element, if maximum wave number selected is less
#' than the wave number length, then maximum wave number selected is coerced to
#' use all wave number unless user defined another index range through \code{index_q_AIUQ}.
#'
#' If \code{model_name} equals 'user_defined', or NA (will coerced to
#' 'user_defined'), then \code{msd_fn} and \code{num_param} need to be provided
#' for parameter estimation.
#'
#' Number of particles \code{M} is set to 50 or automatically extracted from
#' \code{simulation} class for simulated data using \code{simulation} in AIUQ
#' package.
#'
#' By default, using all wave vectors from complete q ring for both \code{AIUQ},
#' unless user defined index range through \code{index_q_AIUQ}.
#'
#' @return Returns an S4 object of class \code{aniso_SAM}.
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
#' # Example 1: Estimation for simulated data
#' set.seed(1)
#' aniso_sim = aniso_simulation(sz=100,len_t=100, model_name="BM",M=100,sigma_bm=c(0.5,0.3))
#' show(aniso_sim)
#' plot_traj(object=aniso_sim)
#' aniso_sam = aniso_SAM(sim_object=aniso_sim, model_name="BM",AIUQ_thr = c(0.999,0))
#' show(aniso_sam)
#' plot_MSD(aniso_sam,msd_truth = aniso_sam@msd_truth)
aniso_SAM <- function(intensity=NA,intensity_str="T_SS_mat",pxsz=1, sz=c(NA,NA),
                      mindt=1,AIUQ_thr=c(1,1),model_name='BM',sigma_0_2_ini=NaN,
                      param_initial=NA,num_optim=1,msd_fn=NA,msd_grad_fn=NA,
                      num_param=NA,uncertainty=FALSE,M=50,sim_object=NA, msd_truth=NA,
                      method="AIUQ",index_q_AIUQ=NA, message_out=TRUE,
                      square=FALSE){

  model <- methods::new("aniso_SAM")
  #check
  if(!is.character(intensity_str)){
    stop("Structure of the intensity profile input should be a character. \n")
  }
  if(intensity_str!="SST_array" && intensity_str!="S_ST_mat" && intensity_str!="T_SS_mat"){
    stop("Structure of the intensity profile input should be one of the type listed in help page. \n")
  }
  if(!is.numeric(pxsz)){
    stop("Pixel size should be a numerical value. \n")
  }
  if(!is.numeric(mindt)){
    stop("Lag time between 2 consecutive image should be a numerical value. \n")
  }

  if(length(AIUQ_thr)==1){
    AIUQ_thr = c(AIUQ_thr,1)
  }
  if(is.na(AIUQ_thr[1])){
    AIUQ_thr = c(1,AIUQ_thr[2])
  }
  if(is.na(AIUQ_thr[2])){
    AIUQ_thr = c(AIUQ_thr[1],1)
  }
  if(!is.numeric(AIUQ_thr)){
    stop("AIUQ threshold should be a numerical vector. \n")
  }
  if(AIUQ_thr[1]<0 || AIUQ_thr[2]<0 || AIUQ_thr[1]>1 || AIUQ_thr[2]>1){
    stop("AIUQ threshold has value between 0 and 1. \n")
  }

  if(!is.numeric(sigma_0_2_ini)){
    stop("Inital value for background noise should be numeric. \n")
  }

  if(class(sim_object)[1]=="simulation"){
    intensity = sim_object@intensity
    model@pxsz = sim_object@pxsz
    model@mindt = sim_object@mindt
    M = sim_object@M
    len_t = sim_object@len_t
    param = matrix(rep(sim_object@param,2),nrow=length(sim_object@param),ncol=2)
    if(sim_object@model_name=="BM"){
      model_param = apply(param,2,function(x){get_true_param_aniso_sim(param_truth=x,model_name=sim_object@model_name)})
      model@param_truth = matrix(model_param,nrow=1,ncol=2)
    }else{
      model@param_truth = apply(param,2,function(x){get_true_param_aniso_sim(param_truth=x,model_name=sim_object@model_name)})
    }
    model@msd_truth = apply(model@param_truth, 2,function(x){get_MSD(theta = x ,d_input=0:(len_t-1),model_name=sim_object@model_name)})
    model@sigma_2_0_truth = sim_object@sigma_2_0
    sz = sim_object@sz
  }else if(class(sim_object)[1]=="aniso_simulation"){
    intensity = sim_object@intensity
    model@pxsz = sim_object@pxsz
    model@mindt = sim_object@mindt
    M = sim_object@M
    len_t = sim_object@len_t
    if(sim_object@model_name=="BM"){
      model_param = apply(sim_object@param,2,function(x){get_true_param_aniso_sim(param_truth=x,model_name=sim_object@model_name)})
      model@param_truth = matrix(model_param,1)
    }else{
      model@param_truth = apply(sim_object@param,2,function(x){get_true_param_aniso_sim(param_truth=x,model_name=sim_object@model_name)})
    }
    model@msd_truth = sim_object@theor_msd
    model@sigma_2_0_truth = sim_object@sigma_2_0
    sz = sim_object@sz
  }else{
    model@pxsz = pxsz
    model@mindt = mindt
    model@msd_truth = matrix(msd_truth,1,1)
    model@sigma_2_0_truth = NA
    model@param_truth = matrix(NA,1,1)
    sz = sz
  }
  if(is.vector(intensity)){
    if(is.na(intensity)){
      stop("Intensity profile can't be missing and should have one of the structure listed in intensity_str. \n")
    }
  }

  if(is.na(model_name)){
    model_name = "user_defined"
  }
  if (model_name == "user_defined" && is.na(num_param)){
    stop("For user defined model, number of parameters that need to be estimated can't be empty. \n")
  }

  if(!is.character(model_name)){
    stop("Fitted model name should be character. \n")
  }
  if(!is.character(method)){
    stop("Method should be character. \n")
  }

  model@model_name = model_name
  model@method = method

  # Transform intensity into the same format and crop image into square image
  # total number of pixels in each image = sz_x*sz_y
  intensity_list = intensity_format_transform(intensity = intensity,
                                              intensity_str = intensity_str,
                                              square = square,sz=sz)

  # Fourier transform
  fft_list = FFT2D(intensity_list=intensity_list,pxsz=model@pxsz,mindt=model@mindt)

  #num of rows and columns of intensity matrix, also representing frame size in y and x directions
  model@sz = c(fft_list$sz_y,fft_list$sz_x)

  model@len_q = fft_list$len_q
  model@len_t = fft_list$len_t
  model@q = fft_list$q
  model@d_input = fft_list$d_input


  if(!is.na(index_q_AIUQ)[1]){
    if(min(index_q_AIUQ)<1 || max(index_q_AIUQ)>model@len_q){
      stop("Selected q range should between 1 and half frame size. \n")
    }
  }



  # get each q ring location index
  if(model@sz[1]==model@sz[2]){
    v = (-(model@sz[1]-1)/2):((model@sz[1]-1)/2)
    x = matrix(rep(v,each = model@sz[1]), byrow = FALSE,nrow = model@sz[1])
    y = matrix(rep(v,each = model@sz[1]), byrow = TRUE,nrow = model@sz[1])
  }else{
    v_x = (-(model@sz[2]-1)/2):((model@sz[2]-1)/2)
    v_y = (-(model@sz[1]-1)/2):((model@sz[1]-1)/2)
    x = matrix(rep(v_x,each = model@sz[1]), byrow = FALSE,nrow = model@sz[1])
    y = matrix(rep(v_y,each = model@sz[2]), byrow = TRUE,nrow = model@sz[1])
  }
  theta_q = cart2polar(x, y)
  q_ring_num = theta_q[,(model@sz[2]+1):dim(theta_q)[2]]
  q_ring_num = round(q_ring_num)


  nq_index = vector(mode = "list")
  for(i in 1:model@len_q){
    nq_index[[i]] = which(q_ring_num==i)
  }
  q_ori_ring_loc = fftshift(q_ring_num, dim = 3)

  q_ori_ring_loc_index = as.list(1:model@len_q)
  total_q_ori_ring_loc_index = NULL
  for(i in 1:model@len_q){
    q_ori_ring_loc_index[[i]] = which(q_ori_ring_loc==i)
    total_q_ori_ring_loc_index = c(total_q_ori_ring_loc_index, q_ori_ring_loc_index[[i]])
  }

  #model@q_ring_loc = q_ring_num
  q_ring_loc = q_ring_num
  #model@q_ori_ring_loc = q_ori_ring_loc
  #model@q_ori_ring_loc_index = q_ori_ring_loc_index
  #model@total_q_ori_ring_loc_index = total_q_ori_ring_loc_index

  #Get A an B ini est
  avg_I_2_ori = 0

  for(i in 1:model@len_t){
    avg_I_2_ori = avg_I_2_ori+abs(fft_list$I_q_matrix[,i])^2/(model@sz[1]*model@sz[2])
  }
  avg_I_2_ori = avg_I_2_ori/model@len_t

  model@I_o_q_2_ori = rep(NA,model@len_q)
  for(i in 1:model@len_q){
    model@I_o_q_2_ori[i] = mean(avg_I_2_ori[q_ori_ring_loc_index[[i]]])
  }
  I_o_q_2_ori_last = model@I_o_q_2_ori[model@len_q]

  model@B_est_ini = 2*I_o_q_2_ori_last
  model@A_est_ini = 2*(model@I_o_q_2_ori - I_o_q_2_ori_last)

  for (i in 1:model@len_q){
    if(sum(model@A_est_ini[1:i])/sum(model@A_est_ini)>=AIUQ_thr[1]){
      num_q_max = i
      break
    }
  }
  if(num_q_max/model@len_q<=AIUQ_thr[2]){
    num_q_max=ceiling(AIUQ_thr[2]*model@len_q)
  }

  # get unique index

  #model@q_ori_ring_loc_unique_index = as.list(1:model@len_q)
  q_ori_ring_loc_unique_index = as.list(1:model@len_q)
  # for(i in 1:model@len_q){
  #   unique_val = unique(avg_I_2_ori[q_ori_ring_loc_index[[i]]])
  #   unique_val = unique_val[1:(length(q_ori_ring_loc_index[[i]])/2)]
  #   index_selected = NULL
  #   for(j in 1:length(unique_val)){
  #     index_selected = c(index_selected,which(avg_I_2_ori == unique_val[j])[1])
  #   }
  #   q_ori_ring_loc_unique_index[[i]] = index_selected
  # }
  for(i in 1:model@len_q){
    len_here = (length(q_ori_ring_loc_index[[i]])-2)/2
    q_ori_ring_loc_unique_index[[i]] = q_ori_ring_loc_index[[i]][c(1,3:(3+len_here-1))]
  }

  total_q_ori_ring_loc_unique_index = NULL
  for(i in 1:model@len_q){
    total_q_ori_ring_loc_unique_index=c(total_q_ori_ring_loc_unique_index,
                                        q_ori_ring_loc_unique_index[[i]])
  }

  #anisotropic
  if(model@sz[1]==model@sz[2]){
    q1_unique_index=as.list(1:model@len_q)
    q2_unique_index=as.list(1:model@len_q)
    for(i in 1:model@len_q){
      index_here=q_ori_ring_loc_unique_index[[i]]
      total_num_unique_index_here=length(q_ori_ring_loc_unique_index[[i]])
      q1_unique_index[[i]]=q2_unique_index[[i]]=rep(NA,total_num_unique_index_here)
      for(j in 1:(total_num_unique_index_here)){
        q1_unique_index[[i]][j]=floor((index_here[j]-1)/model@sz[1]) ##could contain zero
        left_here=(index_here[j]-1)%%model@sz[1]
        if(left_here<=model@len_q){
          q2_unique_index[[i]][j]=(index_here[j]-1)%%model@sz[1] ##could contain zero
        }else{
          q2_unique_index[[i]][j]=model@sz[1]-(index_here[j]-1)%%model@sz[1]-1 ##could contain zero
        }
      }
    }
    q1 = c((1:((model@sz[1]-1)/2))*2*pi/(model@sz[1]*model@pxsz))
    q2 = c((1:((model@sz[1]-1)/2))*2*pi/(model@sz[1]*model@pxsz))
  }else{stop("Update for rectangular image!")}


  if(is.na(sigma_0_2_ini)){sigma_0_2_ini=min(model@I_o_q_2_ori)}

  if(sum(is.na(param_initial))>=1){
    param_initial = get_initial_param(model_name=paste(model@model_name,"_anisotropic",sep=""),
                                      sigma_0_2_ini=sigma_0_2_ini,
                                      num_param=num_param)
  }else{
    param_initial = log(c(param_initial,sigma_0_2_ini))
  }

  if(model@method == "AIUQ"){
    ##if index_q_AIUQ is not defined then we define it; otherwise we use the defined index
    if( is.na(index_q_AIUQ)[1]){
      index_q_AIUQ = 1:num_q_max
    }

    p = (length(param_initial)-1)/2
    num_iteration_max = 50+(2*p-1)*10
    lower_bound = c(rep(-30,2*p),-Inf)
    if(model@model_name == "user_defined"){
      if(is.function(msd_grad_fn)==T){
        gr = anisotropic_log_lik_grad
      }else{gr = NULL}
    }else{gr = anisotropic_log_lik_grad}

    m_param = try(optim(param_initial,anisotropic_log_lik, #gr=gr,
                        I_q_cur=fft_list$I_q_matrix,
                        B_cur=NA,index_q=index_q_AIUQ,
                        I_o_q_2_ori=model@I_o_q_2_ori,
                        q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                        sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                        q1=q1,q2=q2,q1_unique_index=q1_unique_index,
                        q2_unique_index=q2_unique_index,
                        model_name=paste(model@model_name,"_anisotropic",sep=""),
                        msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,
                        method='L-BFGS-B',lower=lower_bound,
                        control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)
    if(num_optim>1){
      for(i_try in 1:(num_optim-1)){
        param_initial_try=param_initial+i_try*runif(2*p+1)
        if(message_out){
          cat("start of another optimization, initial values: ",param_initial_try, "\n")
        }
        m_param_try = try(optim(param_initial_try,anisotropic_log_lik, #gr=gr,
                                I_q_cur=fft_list$I_q_matrix,B_cur=NA,
                                index_q=index_q_AIUQ,I_o_q_2_ori=model@I_o_q_2_ori,
                                q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                                sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                                qq1=q1,q2=q2,q1_unique_index=q1_unique_index,
                                q2_unique_index=q2_unique_index,
                                model_name=paste(model@model_name,"_anisotropic",sep=""),
                                msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,
                                method='L-BFGS-B',lower=lower_bound,
                                control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)
        if(class(m_param)[1]!="try-error"){
          if(class(m_param_try)[1]!="try-error"){
            if(m_param_try$value>m_param$value){
              m_param=m_param_try
            }
          }
        }else{##if m_param has an error then change
          m_param=m_param_try
        }

      }
    }

    count_compute=0 ##if it has an error in optimization, try some more
    while(class(m_param)[1]=="try-error"){
      count_compute=count_compute+1
      param_initial_try=param_initial+count_compute*runif(2*p+1)
      if(message_out){
        cat("start of another optimization, initial values: ",param_initial_try, "\n")
      }
      m_param = try(optim(param_initial_try,anisotropic_log_lik, #gr=gr,
                          I_q_cur=fft_list$I_q_matrix,B_cur=NA,
                          index_q=index_q_AIUQ,I_o_q_2_ori=model@I_o_q_2_ori,
                          q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                          sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                          q1=q1,q2=q2,q1_unique_index=q1_unique_index,
                          q2_unique_index=q2_unique_index,
                          model_name=paste(model@model_name,"_anisotropic",sep=""),
                          msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,method='L-BFGS-B',
                          lower=lower_bound,control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)
      if(count_compute>=2){
        break
      }
    }

    ##if still not converge to a finite value, let's try no derivative search
    if(class(m_param)[1]=="try-error"){
      m_param = try(optim(param_initial,anisotropic_log_lik,
                          I_q_cur=fft_list$I_q_matrix,
                          B_cur=NA,index_q=index_q_AIUQ,
                          I_o_q_2_ori=model@I_o_q_2_ori,
                          q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                          sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                          q1=q1,q2=q2,q1_unique_index=q1_unique_index,
                          q2_unique_index=q2_unique_index,
                          model_name=paste(model@model_name,"_anisotropic",sep=""),
                          msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,method='L-BFGS-B',
                          lower=lower_bound,control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)

      count_compute=0
      while(class(m_param)[1]=="try-error"){
        count_compute=count_compute+1
        #compute_twice=T
        ##change it to runif
        #c(rep(0.5,p),0)
        param_initial_try=param_initial+count_compute*runif(2*p+1)
        if(message_out){
          cat("start of another optimization, initial values: ",param_initial_try, "\n")
        }
        m_param = try(optim(param_initial_try,anisotropic_log_lik,
                            I_q_cur=fft_list$I_q_matrix,B_cur=NA,
                            index_q=index_q_AIUQ,I_o_q_2_ori=model@I_o_q_2_ori,
                            q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                            sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                            q1=q1,q2=q2,q1_unique_index=q1_unique_index,
                            q2_unique_index=q2_unique_index,
                            model_name=paste(model@model_name,"_anisotropic",sep=""),
                            msd_fn=msd_fn,msd_grad_fn=msd_grad_fn,
                            method='L-BFGS-B',lower=lower_bound,
                            control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)
        if(count_compute>=2){
          break
        }
      }
    }

    param_est = m_param$par
    model@mle = m_param$value
    AIC = 2*(length(param_est)+length(index_q_AIUQ)-m_param$value)
    model@sigma_2_0_est = exp(param_est[length(param_est)])

    param_est = matrix(param_est[-length(param_est)],ncol=2)
    if(model_name=="BM"){
      param_est = apply(param_est,2,function(x){get_est_param(theta=exp(x),model_name = model@model_name)})
      model@param_est = matrix(param_est,nrow=1,ncol=2)
    }else{
      model@param_est = apply(param_est,2,function(x){get_est_param(theta=exp(x),model_name = model@model_name)})
    }
    model@msd_est = apply(model@param_est, 2,function(x){get_MSD(theta = x,d_input=model@d_input,model_name=model@model_name, msd_fn=msd_fn)})


    if(uncertainty==TRUE && is.na(M)!=TRUE){
      param_uq_range=param_uncertainty_anisotropic(param_est=m_param$par,I_q_cur=fft_list$I_q_matrix,
                                                   index_q=index_q_AIUQ,
                                                   I_o_q_2_ori=model@I_o_q_2_ori,
                                                   q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                                                   sz=model@sz,len_t=model@len_t,
                                                   q1=q1,q2=q2,q1_unique_index=q1_unique_index,
                                                   q2_unique_index=q2_unique_index,
                                                   d_input=model@d_input,
                                                   model_name=paste(model@model_name,"_anisotropic",sep=""),
                                                   M=M,num_iteration_max=num_iteration_max,
                                                   lower_bound=lower_bound,msd_fn=msd_fn,
                                                   msd_grad_fn=msd_grad_fn)


      for(i_p in 1:length(m_param$par)){

        param_uq_range[,i_p]=c(min((m_param$par[i_p]),param_uq_range[1,i_p]),
                               max((m_param$par[i_p]),param_uq_range[2,i_p]))

      }
      SAM_range_list=get_est_parameters_MSD_SAM_interval_anisotropic(param_uq_range,
                                                                     model_name=model@model_name,
                                                                     d_input=model@d_input, msd_fn=msd_fn)
      model@uncertainty  = uncertainty
      model@msd_x_lower = SAM_range_list$MSD_x_lower
      model@msd_x_upper = SAM_range_list$MSD_x_upper
      model@msd_y_lower = SAM_range_list$MSD_y_lower
      model@msd_y_upper = SAM_range_list$MSD_y_upper
      model@param_uq_range = cbind(SAM_range_list$est_parameters_lower,SAM_range_list$est_parameters_upper)
    }else{
      model@uncertainty = FALSE
      model@msd_x_lower = NA
      model@msd_x_upper = NA
      model@msd_y_lower = NA
      model@msd_y_upper = NA
      model@param_uq_range = matrix(NA,1,1)
    }
  }

  if(model@method=='AIUQ'){
    model@index_q = index_q_AIUQ
  }

  model@I_q = fft_list$I_q_matrix
  #model@p = p
  model@AIC = AIC
  model@q_ori_ring_loc_unique_index = q_ori_ring_loc_unique_index
  return(model)
}
