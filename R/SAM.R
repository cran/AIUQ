#' Scattering analysis of microscopy
#'
#' @description
#' Fast parameter estimation in scattering analysis of microscopy, using either
#' AIUQ or DDM method.
#'
#' @param intensity intensity profile. See 'Details'.
#' @param intensity_str structure of the intensity profile, options from
#' ('SST_array','S_ST_mat','T_SS_mat'). See 'Details'.
#' @param pxsz size of one pixel in unit of micron, 1 for simulated data
#' @param sz frame size of the intensity profile in x and y directions,
#' number of pixels contained in each frame equals sz_x by sz_y.
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
#' @param method methods for parameter estimation, options from ('AIUQ','DDM_fixedAB','DDM_estAB').
#' @param index_q_AIUQ index range for wave number when using AIUQ method. See 'Details'.
#' @param index_q_DDM index range for wave number when using DDM method. See 'Details'.
#' @param message_out a logical evaluating to TRUE or FALSE indicating whether
#' or not to output the message.
#' @param A_neg controls modification for negative A(q), options from  ('abs','zero'),
#' with setting negative A(q) to its absolute value as the default.
#' @param square a logical evaluating to TRUE or FALSE indicating whether or not
#' to crop the original intensity profile into square image.
#' @param output_dqt a logical evaluating to TRUE or FALSE indicating whether or
#' not to compute observed dynamic image structure function(Dqt).
#' @param output_isf a logical evaluating to TRUE or FALSE indicating whether or
#' not to compute empirical intermediate scattering function(ISF).
#' @param output_modeled_isf a logical evaluating to TRUE or FALSE indicating
#' whether or not to compute modeled intermediate scattering function(ISF).
#' @param output_modeled_dqt a logical evaluating to TRUE or FALSE indicating
#' whether or not to compute modeled dynamic image structure function(Dqt).
#'
#' @details
#' For simulated data using \code{simulation} in AIUQ package, \code{intensity}
#' will be automatically extracted from \code{simulation} class.
#'
#' By default \code{intensity_str} is set to 'T_SS_mat', a time by space\eqn{\times}{%\times}space
#' matrix, which is the structure of intensity profile obtained from \code{simulation}
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
#' By default, using all wave vectors from complete q ring, unless user defined
#' index range through \code{index_q_AIUQ} or \code{index_q_DDM}.
#'
#' @return Returns an S4 object of class \code{SAM}.
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
#' sim_bm = simulation(len_t=100,sz=100,sigma_bm=0.5)
#' show(sim_bm)
#' sam = SAM(sim_object = sim_bm)
#' show(sam)
SAM <- function(intensity=NA,intensity_str="T_SS_mat",pxsz=1,sz=c(NA,NA),mindt=1,
                AIUQ_thr=c(1,1),model_name='BM',sigma_0_2_ini=NaN,param_initial=NA,
                num_optim=1,msd_fn=NA,msd_grad_fn=NA,num_param=NA,
                uncertainty=FALSE,M=50,sim_object=NA, msd_truth=NA,
                method="AIUQ",index_q_AIUQ=NA, index_q_DDM=NA,message_out=TRUE,
                A_neg="abs", square=FALSE, output_dqt=FALSE, output_isf=FALSE,
                output_modeled_isf=FALSE,output_modeled_dqt=FALSE){

  model <- methods::new("SAM")
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
    #model@mindt = sim_object@mindt
    M = sim_object@M
    #model_name = sim_object@model_name ##not always inherent
    len_t = sim_object@len_t
    model@param_truth = get_true_param_sim(param_truth=sim_object@param,model_name=sim_object@model_name)
    model@msd_truth = sim_object@theor_msd
    model@sigma_2_0_truth = sim_object@sigma_2_0
    sz = sim_object@sz
  }else{
    model@pxsz = pxsz
    #model@mindt = mindt
    model@msd_truth = msd_truth
    model@sigma_2_0_truth = NA
    model@param_truth = NA
    sz = sz
  }
  model@mindt = mindt
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

  if(!is.na(index_q_DDM)[1]){
    if(min(index_q_DDM)<1 || max(index_q_DDM)>model@len_q){
      stop("Selected q range should between 1 and half frame size. \n")
    }
  }

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

  B_est_ini = 2*I_o_q_2_ori_last
  A_est_ini = 2*(model@I_o_q_2_ori - I_o_q_2_ori_last)

  for (i in 1:model@len_q){
    if(sum(A_est_ini[1:i])/sum(A_est_ini)>=AIUQ_thr[1]){
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


  if(is.na(sigma_0_2_ini)){sigma_0_2_ini=min(model@I_o_q_2_ori)}
  if(sum(is.na(param_initial))>=1){
    param_initial = get_initial_param(model_name=model@model_name,
                                      sigma_0_2_ini=sigma_0_2_ini,
                                      num_param=num_param)
  }else{
    param_initial = c(param_initial,sigma_0_2_ini)
  }

  if(model@method == "AIUQ"){
    ##if index_q_AIUQ is not defined then we define it; otherwise we use the defined index
    if( is.na(index_q_AIUQ)[1]){
      index_q_AIUQ = 1:num_q_max
    }

    p = length(param_initial)-1
    num_iteration_max = 50+(p-1)*10
    lower_bound = c(rep(-30,p),-Inf)
    if(model@model_name == "user_defined"){
      if(is.function(msd_grad_fn)==T){
        gr = log_lik_grad
      }else{gr = NULL}
    }else if(A_neg=="zero"){
      gr = NULL
    }else{gr = log_lik_grad}
    m_param = try(optim(param_initial,log_lik, gr=gr,
                        I_q_cur=fft_list$I_q_matrix,
                        B_cur=NA, A_neg = A_neg, index_q=index_q_AIUQ,
                        I_o_q_2_ori=model@I_o_q_2_ori,
                        q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                        sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                        q=model@q,model_name=model@model_name,msd_fn=msd_fn,
                        msd_grad_fn=msd_grad_fn,method='L-BFGS-B',lower=lower_bound,
                        control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)
     if(num_optim>1){
       for(i_try in 1:(num_optim-1)){
         param_initial_try=param_initial+i_try*runif(p+1)
         if(message_out){
           cat("start of another optimization, initial values: ",param_initial_try, "\n")
         }
         m_param_try = try(optim(param_initial_try,log_lik, gr=gr,
                             I_q_cur=fft_list$I_q_matrix,
                             B_cur=NA, A_neg = A_neg,index_q=index_q_AIUQ,
                             I_o_q_2_ori=model@I_o_q_2_ori,
                             q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                             sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                             q=model@q,model_name=model@model_name,msd_fn=msd_fn,
                             msd_grad_fn=msd_grad_fn,method='L-BFGS-B',lower=lower_bound,
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
      param_initial_try=param_initial+count_compute*runif(p+1)
      if(message_out){
        cat("start of another optimization, initial values: ",param_initial_try, "\n")
      }
      m_param = try(optim(param_initial_try,log_lik,gr=gr,
                          I_q_cur=fft_list$I_q_matrix,B_cur=NA, A_neg = A_neg,
                          index_q=index_q_AIUQ,I_o_q_2_ori=model@I_o_q_2_ori,
                          q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                          sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                          q=model@q,model_name=model@model_name,msd_fn=msd_fn,
                          msd_grad_fn=msd_grad_fn,method='L-BFGS-B',lower=lower_bound,
                          control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)
      if(count_compute>=2){
        break
      }
    }

    ##if still not converge to a finite value, let's try no derivative search
    if(class(m_param)[1]=="try-error"){
      m_param = try(optim(param_initial,log_lik,
                          I_q_cur=fft_list$I_q_matrix,
                          B_cur=NA, A_neg = A_neg,index_q=index_q_AIUQ,
                          I_o_q_2_ori=model@I_o_q_2_ori,
                          q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                          sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                          q=model@q,model_name=model@model_name,msd_fn=msd_fn,
                          msd_grad_fn=msd_grad_fn,method='L-BFGS-B',lower=lower_bound,
                          control = list(fnscale=-1,maxit=num_iteration_max)),TRUE)

      count_compute=0
      while(class(m_param)[1]=="try-error"){
        count_compute=count_compute+1
        #compute_twice=T
        ##change it to runif
        #c(rep(0.5,p),0)
        param_initial_try=param_initial+count_compute*runif(p+1)
        if(message_out){
          cat("start of another optimization, initial values: ",param_initial_try, "\n")
        }
        m_param = try(optim(param_initial_try,log_lik,
                            I_q_cur=fft_list$I_q_matrix,B_cur=NA, A_neg = A_neg,
                            index_q=index_q_AIUQ,I_o_q_2_ori=model@I_o_q_2_ori,
                            q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                            sz=model@sz,len_t=model@len_t,d_input=model@d_input,
                            q=model@q,model_name=model@model_name,msd_fn=msd_fn,
                            msd_grad_fn=msd_grad_fn,method='L-BFGS-B',lower=lower_bound,
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
    #est_list = get_est_param_MSD(theta=exp(param_est),d_input=model@d_input,model_name = model@model_name)

    model@param_est = get_est_param(theta=exp(param_est),model_name = model@model_name)
    if(model@model_name=='user_defined'){
      model@param_est = model@param_est[-length(param_initial)]
    }
    model@msd_est = get_MSD(theta=model@param_est,d_input=model@d_input,
                            model_name = model@model_name, msd_fn=msd_fn)

    if(uncertainty==T && is.na(M)!=T){
      param_uq_range=param_uncertainty(param_est=m_param$par,I_q_cur=fft_list$I_q_matrix,
                                       index_q=index_q_AIUQ,A_neg=A_neg,
                                       I_o_q_2_ori=model@I_o_q_2_ori,
                                       q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                                       sz=model@sz,len_t=model@len_t,q=model@q,
                                       d_input=model@d_input,
                                       model_name=model@model_name,M=M,
                                       num_iteration_max=num_iteration_max,
                                       lower_bound=lower_bound,msd_fn=msd_fn,
                                       msd_grad_fn=msd_grad_fn)


      for(i_p in 1:length(m_param$par)){

        param_uq_range[,i_p]=c(min((m_param$par[i_p]),param_uq_range[1,i_p]),
                               max((m_param$par[i_p]),param_uq_range[2,i_p]))

      }
      SAM_range_list=get_est_parameters_MSD_SAM_interval(param_uq_range,
                                                         model_name=model@model_name,
                                                         d_input=model@d_input, msd_fn=msd_fn)
      model@uncertainty  = uncertainty
      model@msd_lower = SAM_range_list$MSD_lower
      model@msd_upper = SAM_range_list$MSD_upper
      model@param_uq_range = cbind(SAM_range_list$est_parameters_lower,SAM_range_list$est_parameters_upper)
    }else{
      model@uncertainty  = uncertainty
      model@msd_lower = NA
      model@msd_upper = NA
      model@param_uq_range = matrix(NA,1,1)
    }

    if(output_dqt==FALSE && output_isf==FALSE){
      Dqt = matrix(NA,1,1)
      isf = matrix(NA,1,1)
    }else if(output_dqt==TRUE && output_isf==FALSE){
      Dqt = SAM_Dqt(len_q=model@len_q,index_q=1:model@len_q,len_t=model@len_t,
                    I_q_matrix=fft_list$I_q_matrix,sz=model@sz,
                    q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index)
      isf = matrix(NA,1,1)
    }else if(output_isf==TRUE){
      Dqt = matrix(NA,model@len_q,model@len_t-1)
      isf = matrix(NA,model@len_q,model@len_t-1)
      for (q_j in 1:model@len_q){
        index_cur = q_ori_ring_loc_unique_index[[q_j]]
        I_q_cur = fft_list$I_q_matrix[index_cur,]
        for (t_i in 1:(model@len_t-1)){
          Dqt[q_j,t_i]=mean((abs(I_q_cur[,(t_i+1):model@len_t]-I_q_cur[,1:(model@len_t-t_i)]))^2/(model@sz[1]*model@sz[2]),na.rm=T)
        }
        if(A_est_ini[q_j]==0){break}
        isf[q_j,] = 1-(Dqt[q_j,]-B_est_ini)/A_est_ini[q_j]
      }
    }

    model@B_est = model@sigma_2_0_est*2
    if(A_neg=="abs"){
      model@A_est = abs(2*(model@I_o_q_2_ori - model@B_est/2))
    }else if(A_neg=="zero"){
      model@A_est = 2*(model@I_o_q_2_ori- model@B_est/2)
      model@A_est = ifelse(model@A_est>0,model@A_est,0)
    }

  }else if(model@method == "DDM_fixedAB"){
    if(length(index_q_DDM)==1 &&is.na(index_q_DDM)){
      index_q_DDM = 1:model@len_q
    }
    Dqt = SAM_Dqt(len_q=model@len_q,index_q=index_q_DDM,len_t=model@len_t,
                  I_q_matrix=fft_list$I_q_matrix,sz=model@sz,
                  q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index)
    if(output_isf==TRUE){
      isf = matrix(NA,model@len_q,model@len_t-1)
      for (q_j in 1:model@len_q){
        if(A_est_ini[q_j]==0){break}
        isf[q_j,] = 1-(Dqt[q_j,]-B_est_ini)/A_est_ini[q_j]
        }
      }else{
        isf = matrix(NA,1,1)
      }

    l2_est_list = theta_est_l2_dqt_fixedAB(param=param_initial[-length(param_initial)],q=model@q,index_q=index_q_DDM,
                                           Dqt=Dqt,A_est_q=A_est_ini,B_est=B_est_ini,
                                           d_input=model@d_input, model_name=model@model_name,
                                           msd_fn=msd_fn,msd_grad_fn=msd_grad_fn)


    model@param_est = l2_est_list$param_est
    model@msd_est = l2_est_list$msd_est
    model@sigma_2_0_est = B_est_ini/2
    model@A_est = l2_est_list$A_est
    p = NaN
    AIC = NaN
    model@mle = NaN
    model@param_uq_range = matrix(NA,1,1)
  }else if(model@method == "DDM_estAB"){
    if(length(index_q_DDM)==1 &&is.na(index_q_DDM)){
      index_q_DDM = 1:model@len_q
    }
    Dqt = SAM_Dqt(len_q=model@len_q,index_q=index_q_DDM,len_t=model@len_t,
                  I_q_matrix=fft_list$I_q_matrix,sz=model@sz,
                  q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index)
    if(output_isf==TRUE){
      isf = matrix(NA,model@len_q,model@len_t-1)
      for (q_j in 1:model@len_q){
        if(A_est_ini[q_j]==0){break}
        isf[q_j,] = 1-(Dqt[q_j,]-B_est_ini)/A_est_ini[q_j]
      }
    }else{
      isf = matrix(NA,1,1)
    }

    l2_est_list = theta_est_l2_dqt_estAB(param=param_initial,q=model@q,index_q=index_q_DDM,
                                         Dqt=Dqt,A_ini=A_est_ini,d_input=model@d_input,
                                         model_name=model@model_name,msd_fn=msd_fn,msd_grad_fn=msd_grad_fn)

    model@param_est = l2_est_list$param_est
    model@msd_est = l2_est_list$msd_est
    model@sigma_2_0_est = l2_est_list$sigma_2_0_est
    model@A_est = l2_est_list$A_est
    p = NaN
    AIC = NaN
    model@mle = NaN
    model@param_uq_range = matrix(NA,1,1)
  }

  model@Dqt = Dqt
  model@ISF = isf

  if(output_modeled_dqt==FALSE && output_modeled_isf==FALSE){
    model@modeled_Dqt = matrix(NA,1,1)
    model@modeled_ISF = matrix(NA,1,1)
  }else if(output_modeled_isf==TRUE && output_modeled_dqt==FALSE){
    model@modeled_ISF = matrix(NA,model@len_q,model@len_t-1)
    model@modeled_Dqt = matrix(NA,1,1)
    for(q_j in 1:model@len_q){
      q_selected = model@q[q_j]
      model@modeled_ISF [q_j,] = exp(-q_selected^2*model@msd_est[-1]/4)
    }
  }else if(output_modeled_dqt==TRUE){
    model@modeled_ISF  = matrix(NA,model@len_q,model@len_t-1)
    model@modeled_Dqt = matrix(NA,model@len_q,model@len_t-1)
    for(q_j in 1:model@len_q){
      q_selected = model@q[q_j]
      model@modeled_ISF[q_j,] = exp(-q_selected^2*model@msd_est[-1]/4)
      if(A_est_ini[q_j]==0){break}
      model@modeled_Dqt[q_j,] = A_est_ini[q_j]*(1-model@modeled_ISF[q_j,])+model@sigma_2_0_est*2
    }

  }

  if(model@method=='AIUQ'){
     model@index_q = index_q_AIUQ
  }else{ ##DDM
    model@index_q = index_q_DDM

  }
  model@I_q = fft_list$I_q_matrix
  #model@p = p
  model@AIC = AIC
  model@q_ori_ring_loc_unique_index = q_ori_ring_loc_unique_index
  return(model)
}
