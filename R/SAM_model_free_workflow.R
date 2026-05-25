
SAM_model_free_workflow <- function(model, intensity = NA, sim_object = NA, intensity_str = NA,
                                    mindt = 1, pxsz = 1, uncertainty = FALSE, M = NA,
                                    estimate_moduli = FALSE,
                                    MSD_unit = 1, kB = 1.380649*10^-23,
                                    Temp = NA, a = NA, substack = FALSE, 
                                    ...) {
  
  
  if (class(sim_object)[1]=="simulation") {
    # simulation
    rng = range(sim_object@intensity, na.rm = TRUE)
    if (isTRUE(all.equal(rng[1], 0)) && isTRUE(all.equal(rng[2], 1))) {
      raw_intensity = sim_object@intensity
    } else {
      raw_intensity = (sim_object@intensity - rng[1]) / (rng[2] - rng[1])
    }
    
    dim_intensity = c(sim_object@sz[1], sim_object@sz[2], sim_object@len_t)
    current_format = "T_SS_mat" 
    M = sim_object@M
    pxsz = sim_object@pxsz
    mindt = sim_object@mindt
    model@msd_truth = sim_object@theor_msd
    #model@sigma_2_0_truth = sim_object@sigma_2_0
  } else {
    # real data
    rng = range(intensity, na.rm = TRUE)
    if (isTRUE(all.equal(rng[1], 0)) && isTRUE(all.equal(rng[2], 1))) {
      raw_intensity = intensity
    } else {
      raw_intensity = (intensity - rng[1]) / (rng[2] - rng[1])
    }
    
    rm(intensity)
    gc()
    
    dim_intensity = dim(raw_intensity)
    if (is.null(dim_intensity) || length(dim_intensity) != 3) {
      stop("Input intensity should have size space x space x time. \n")
    }
    current_format = intensity_str
    model@msd_truth = NA
    #model@sigma_2_0_truth = NA
  }
  
  q_thr = 0.999
  data_fft = FFT2D_mf(intensity = raw_intensity, ori_format = current_format, pxsz = pxsz, mindt = mindt)
  rm(raw_intensity)
  gc()
  
  
  index_list=Get_q_ring_loc_mf(data_fft$sz, data_fft$len_q)
  ini_est_list=Get_A_B_ini_num_q_max_est_mf(sz = data_fft$sz, len_q = data_fft$len_q, len_t = data_fft$len_t,
                                         I_q_matrix = data_fft$I_q_matrix,
                                         q_ori_ring_loc_index = index_list$q_ori_ring_loc_index,
                                         threshold = q_thr, beta = 0) 
  
  
  num_selected_upper = 20
  upper_limits = log(ini_est_list$num_q_max - ini_est_list$num_q_max/(2*num_selected_upper))
  q_index_selected_here = unique(ceiling(exp(seq(from = log(1), to = upper_limits, 
                                                 by = upper_limits/(num_selected_upper-1))))) 
  
  num_iteration_max = 100
  subsample_t = floor(data_fft$len_t/100)
  
  param_list = list()
  param_list$mindt = mindt
  param_list$pxsz = pxsz
  param_list$nframes = dim_intensity[3]
  param_list$nx = dim_intensity[1]
  param_list$ny = dim_intensity[2]
  param_list$len_q = data_fft$len_q
  param_list$q = data_fft$q
  param_list$I_q_matrix = data_fft$I_q_matrix
  param_list$I_o_q_2 = ini_est_list$I_o_q_2_ori
  param_list$q_ori_ring_loc_index = index_list$q_ori_ring_loc_index
  
  DDM_initial_MSD = Get_initial_param_DDM(param_list)
  
  num_estim = 6
  slope_factor = 0.9

  design_optimization = 'log_equal_space' 
  model_name = 'direct_nonparametric'
  initial_param = get_initial_param_nonparametric_mf(model_name = model_name, d_input = data_fft$d_input, 
                                                  num_estim = num_estim, sigma_0_2_ini = min(ini_est_list$I_o_q_2_ori),
                                                  design_optimization = design_optimization)
  
  
  param_ini_DDM = vector(length = num_estim+1)
  param_ini_DDM[1] = log(DDM_initial_MSD$MSD[1]) 
  param_ini_DDM[length(param_ini_DDM)] = initial_param$param_initial[num_estim+1]
  for (i in 1:num_estim){
    param_ini_DDM[i] = param_ini_DDM[1] + slope_factor * (initial_param$input_training[i]-log(mindt))*((log(DDM_initial_MSD$MSD[2])-log(DDM_initial_MSD$MSD[1]))/(log(2*mindt)-log(mindt)))
  }
  
  # optimization ----
    # stage 1
  loglik_fn = Get_loglik_nonparametric_mf(I_q_cur = data_fft$I_q_matrix, B_cur = NA, num_q_cur = ini_est_list$num_q_max,
                                          I_o_q_2_ori = ini_est_list$I_o_q_2_ori, q_ori_ring_loc_unique_index = ini_est_list$q_ori_ring_loc_unique_index,
                                          sz = data_fft$sz, len_t = data_fft$len_t, d_input = data_fft$d_input, q = data_fft$q,
                                          model_name = "direct_nonparametric", input_training = initial_param$input_training,
                                          q_index_selected = q_index_selected_here, subsample_t = subsample_t
  )
    
  m_param_optim_no_restart_direct_DDM = optim(
    param_ini_DDM, loglik_fn, control = list(fnscale = -1, maxit = 100), method = "L-BFGS-B"
  )
    
  # stage 2
  loglik_fn_fine = Get_loglik_nonparametric_mf(I_q_cur = data_fft$I_q_matrix, B_cur = NA, num_q_cur = ini_est_list$num_q_max,
                                               I_o_q_2_ori = ini_est_list$I_o_q_2_ori, q_ori_ring_loc_unique_index = ini_est_list$q_ori_ring_loc_unique_index,
                                               sz = data_fft$sz, len_t = data_fft$len_t, d_input = data_fft$d_input, q = data_fft$q,
                                               model_name = "direct_nonparametric", input_training = initial_param$input_training,
                                               q_index_selected = q_index_selected_here, substack = substack
  )
    
  m_param_optim_no_restart_direct_DDM_fine = optim(
    m_param_optim_no_restart_direct_DDM$par, loglik_fn_fine, control = list(fnscale = -1, maxit = 50), method = "L-BFGS-B"
  )
  theta_est_optim_no_restart_direct_DDM_fine = exp(m_param_optim_no_restart_direct_DDM_fine$par[-length(m_param_optim_no_restart_direct_DDM_fine$par)])
  MSD_est_optim_no_restart_direct_DDM_fine = Get_MSD_nonparametric_mf(theta = theta_est_optim_no_restart_direct_DDM_fine, d_input = data_fft$d_input, input_training = initial_param$input_training, model_name = "direct_nonparametric")  
  
  
  model@msd_est = MSD_est_optim_no_restart_direct_DDM_fine
  
  if (isTRUE(uncertainty) && !is.na(M)) {
    uq_res=param_uncertainty_nonparametric_mf(param_est=m_param_optim_no_restart_direct_DDM_fine$par,I_q_cur=data_fft$I_q_matrix,B_cur=NA,
                                           num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
                                           q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
                                           sz=data_fft$sz,len_t=data_fft$len_t,
                                           d_input=data_fft$d_input,q=data_fft$q,model_name=model_name,input_training=initial_param$input_training,M=M,
                                           q_index_selected=q_index_selected_here,
                                           num_iteration_max=num_iteration_max, estimation_method = "asymptotics") 
    model@msd_lower = uq_res$MSD_lower
    model@msd_upper = uq_res$MSD_upper
  } else {
    model@msd_lower = NA
    model@msd_upper = NA
  }
  model@d_input = data_fft$d_input
  model@len_t = data_fft$len_t
  model@len_q = data_fft$len_q
  model@index_q = q_index_selected_here
  #model@sigma_2_0_est = exp(m_param_optim_no_restart_direct_DDM_fine$par[length(m_param_optim_no_restart_direct_DDM_fine$par)])
  model@mle = m_param_optim_no_restart_direct_DDM_fine$value
  model@pxsz = pxsz
  model@mindt = mindt
  
  if (isTRUE(estimate_moduli) && !isTRUE(uncertainty)) {
    stop("Modulus estimation requires uncertainty = TRUE because MSD_lower and MSD_upper are needed. \n")
  }
  
  if (isTRUE(estimate_moduli) && isTRUE(uncertainty)) {
    if (!methods::is(sim_object, "simulation")) {
      if (is.na(MSD_unit) || is.na(Temp) || is.na(a)) {
        stop("For modulus estimation using real data, MSD_unit, Temp, and a must be provided. \n")
      }
    }
    
    moduli_res = Get_Gp_Gpp_GSER(MSD_est = MSD_est_optim_no_restart_direct_DDM_fine, d_input = data_fft$d_input,
                                 MSD_lower = model@msd_lower, MSD_upper = model@msd_upper, sim_object = sim_object,
                                 MSD_unit = MSD_unit, kB = kB, Temp = Temp, a = a)
    
    model@omega = moduli_res$omega
    model@Gp = moduli_res$Gprime_est
    model@Gpp = moduli_res$Gdoubleprime_est
    
    if (methods::is(sim_object, "simulation")) {
      model@theor_Gp = moduli_res$Gprime_truth
      model@theor_Gpp = moduli_res$Gdoubleprime_truth
    } else {
      model@theor_Gp = NA
      model@theor_Gpp = NA
    }
  } else {
    model@omega = NA
    model@theor_Gp = NA
    model@theor_Gpp = NA
    model@Gp = NA
    model@Gpp = NA
  }
  
  model@method = "model-free"
  
  return(model)
}






