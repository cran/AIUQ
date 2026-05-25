

intensity_format_transform_mf<-function(ori_format,intensity,edge=0,exclude=0){
  if(ori_format=='SST_array'){ ##two space indices and one time indice, tiff
    #sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    if(dim(intensity)[1]%%2==0){ ##is it correct if it is a rectangle
      sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    }else{
      sz = dim(intensity)[1]
    }
    len_t = dim(intensity)[3] #relevant to format
    
    intensity_transform=matrix(NA,sz^2,len_t)
    for(i in 1:len_t){
      intensity_transform[,i] = as.vector(intensity[1:sz,1:sz,i])  # 1:479; 481:480*2-1,...
    }
    #mid_size = dim(intensity)[1]/2 ##is it useful? 
  }else if(ori_format=='S_ST_mat'){ ##time,  space and time as a matrix, matlab simulation 
    intensity=as.matrix(intensity)
    
    #sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    if(dim(intensity)[1]%%2==0){
      sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    }else{
      sz = dim(intensity)[1]
    }
    
    len_t = dim(intensity)[2]/dim(intensity)[1] #relevant to format
    intensity_transform=matrix(NA,sz^2,len_t)
    for(i in 1:len_t){
      intensity_transform[,i]= as.vector(intensity[1:sz,(1+(sz+1)*(i-1)):(sz+(sz+1)*(i-1))]) 
    }
    
    
  }else if(ori_format=='T_SS_mat'){ ##simulation for R, square matrix, need to think about 
    intensity=as.matrix(intensity)
    
    if(sqrt(dim(intensity)[2])%%2==0){
      sz = sqrt(dim(intensity)[2])-1 #pixel dimension; total # of pixels in each image = sz^2
    }else{
      sz = sqrt(dim(intensity)[2])
    }
    len_t = dim(intensity)[1]
    
    #I_q = array(NaN,dim = c(sz,sz,len_t)) #Fourier transformed intensity in t, not dt
    intensity_transform=matrix(NA,sz^2,len_t)
    for(i in 1:len_t){
      intensity_mat=matrix(intensity[i,],sz+1,sz+1 ) 
      intensity_transform[,i]= as.vector(intensity_mat[1:sz,1:sz] )     
    }
  }
  return(intensity_transform)
}

##having really use edge and exclude 
# FFT2D_mf<-function(intensity,pxsz,mindt){
#   sz=sqrt(dim(intensity)[1])
#   len_t=dim(intensity)[2]
#   I_q_matrix=matrix(NA,sz^2,len_t)
#   for(i in 1:len_t){
#     I_q_matrix[,i] = as.vector(fftwtools::fftw2d(matrix(intensity[,i],sz,sz)))  #
#   }
#   
#   
#   ans_list=list()
#   ans_list$sz=sz
#   ans_list$len_q=length(1:((sz-1)/2))
#   ans_list$len_t=len_t
#   ans_list$I_q_matrix=I_q_matrix
#   ans_list$q= (1:((sz-1)/2))*2*pi/(sz*pxsz)
#   ans_list$input= mindt*(1:(len_t))  ##t
#   ans_list$d_input = ans_list$input[1:length(ans_list$input)]-ans_list$input[1] ##delta t, including zero 
#   
#   return(ans_list)
# }

FFT2D_mf <- function(intensity, ori_format, pxsz, mindt) {
  if (ori_format == "SST_array") {
    if (dim(intensity)[1] %% 2 == 0) {
      sz = dim(intensity)[1] - 1
    } else {
      sz = dim(intensity)[1]
    }
    len_t = dim(intensity)[3]
    
    get_frame <- function(i) {
      intensity[1:sz, 1:sz, i]
    }
    
  } else if (ori_format == "S_ST_mat") {
    intensity = as.matrix(intensity)
    
    if (dim(intensity)[1] %% 2 == 0) {
      sz = dim(intensity)[1] - 1
    } else {
      sz = dim(intensity)[1]
    }
    
    len_t = dim(intensity)[2] / dim(intensity)[1]
    
    get_frame <- function(i) {
      intensity[1:sz, (1 + (sz + 1) * (i - 1)):(sz + (sz + 1) * (i - 1))]
    }
    
  } else if (ori_format == "T_SS_mat") {
    intensity = as.matrix(intensity)
    
    if (sqrt(dim(intensity)[2]) %% 2 == 0) {
      sz = sqrt(dim(intensity)[2]) - 1
    } else {
      sz = sqrt(dim(intensity)[2])
    }
    
    len_t = dim(intensity)[1]
    
    get_frame <- function(i) {
      intensity_mat = matrix(intensity[i, ], sz + 1, sz + 1)
      intensity_mat[1:sz, 1:sz]
    }
    
  } else {
    stop("Unsupported ori_format")
  }
  
  I_q_matrix = matrix(complex(real = 0, imaginary = 0), sz^2, len_t)
  
  for (i in 1:len_t) {
    frame_mat = get_frame(i)
    I_q_matrix[, i] = as.vector(fftwtools::fftw2d(frame_mat))
  }
  
  ans_list = list()
  ans_list$sz = sz
  ans_list$len_q = length(1:((sz - 1) / 2))
  ans_list$len_t = len_t
  ans_list$I_q_matrix = I_q_matrix
  ans_list$q = (1:((sz - 1) / 2)) * 2 * pi / (sz * pxsz)
  ans_list$input = mindt * (1:len_t)
  ans_list$d_input = ans_list$input - ans_list$input[1]
  
  ans_list
}

## Define fftshift: 
# Function that used to swap the first quadrant with the third, and the second quadrant
# with the forth. So the zero-frequency component will be at the center of the array.
fftshift_mf <- function(input_matrix, dim = -1) {
  
  rows <- dim(input_matrix)[1]    
  cols <- dim(input_matrix)[2]    
  
  swap_up_down <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
  }
  
  swap_left_right <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
  }
  
  swap_up_down_reverse <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    return(rbind(input_matrix[((rows_half):rows), (1:cols)], input_matrix[1:(rows_half-1), (1:cols)]))
  }
  
  swap_left_right_reverse <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half):cols)], input_matrix[1:rows, 1:(cols_half-1)]))
  }
  
  
  if (dim == -1) {
    input_matrix <- swap_up_down(input_matrix)
    return(swap_left_right(input_matrix))
  }
  else if (dim == 1) {
    return(swap_up_down(input_matrix))
  }
  else if (dim == 2) {
    return(swap_left_right(input_matrix))
  }else if(dim == 3){
    input_matrix <- swap_up_down_reverse(input_matrix)
    return(swap_left_right_reverse(input_matrix))
    
  }
  else {
    stop("Invalid dimension parameter")
  }
}

Get_q_ring_loc_mf<-function(sz,len_q){
  v = (-(sz-1)/2):((sz-1)/2)
  x = matrix(rep(v,each = sz), byrow = FALSE,nrow=sz)
  y = matrix(rep(v,each = sz), byrow = TRUE,nrow=sz)
  
  theta_q = geometry::cart2pol(x, y)
  #dim(theta_q)
  q_ring_num = theta_q[,(sz+1):dim(theta_q)[2]]
  q_ring_num = round(q_ring_num) ###after transformed q number
  #len_q = length(1:((sz-1)/2))
  nq_index = vector(mode = "list")
  for(i in 1:len_q){
    nq_index[[i]] = which(q_ring_num==i)
  }
  #q_ring_num[nq_index[[3]]] #check whether index is correct 
  #> which(q_ring_num==0)
  #[1] 4901
  
  q_ori_ring_loc=fftshift_mf(q_ring_num, dim = 3)  ###this is the location for original 
  #> which(q_ori_ring_loc==0)
  #[1] 1
  #image2D(q_ori_ring_loc)
  ###only use 1 zero
  # if(length(which(q_ring_num==0))>1){
  #   index_0_num=floor(length(q_ring_num)/2)+1
  #   q_ring_num_0_index=which(q_ring_num==0)
  #   selected_0_index= which(q_ring_num_0_index==index_0_num)
  #   q_ring_num[q_ring_num_0_index[-selected_0_index]]=1 ###change to 1
  #   
  #   q_ori_ring_loc_index=which(q_ori_ring_loc==0)
  #   q_ori_ring_loc[q_ori_ring_loc_index[-1]]=1
  # }
  
  #max(q_ring_num)
  q_ori_ring_loc_index=as.list(1:len_q)
  total_q_ori_ring_loc_index=NULL
  for(i in 1:len_q){
    q_ori_ring_loc_index[[i]]=which(q_ori_ring_loc==i)
    total_q_ori_ring_loc_index=c(total_q_ori_ring_loc_index,  q_ori_ring_loc_index[[i]])
  }
  
  
  
  #image2D(q_ring_num)
  #image2D(q_ori_ring_loc)
  ans_list=list()
  ans_list$q_ring_loc=q_ring_num
  ans_list$q_ori_ring_loc=q_ori_ring_loc
  ans_list$q_ori_ring_loc_index=q_ori_ring_loc_index ##original index 
  ans_list$total_q_ori_ring_loc_index=total_q_ori_ring_loc_index ##original index 
  
  
  return(ans_list)
}


##q_ori_ring_loc_index already contains isotropic information 
Get_A_B_ini_num_q_max_est_mf<-function(sz,len_q,len_t, 
                                    I_q_matrix,q_ori_ring_loc_index,threshold=0.999,beta=0.5,anisotropic=T){
  
  #avg_I_2_ori=matrix(0,nrow = sz, ncol = sz)
  avg_I_2_ori=0
  
  for(i in 1:len_t){
    #avg_I_2_ori=avg_I_2_ori+abs(I_q[,,i])^2/sz^2
    avg_I_2_ori=avg_I_2_ori+abs(I_q_matrix[,i])^2/sz^2
    
  }
  avg_I_2_ori = avg_I_2_ori/len_t
  
  # I_o_q_2_ori=rep(NA,len_q)
  # for(i in 1:len_q){
  #   I_o_q_2_ori[i]=mean(avg_I_2_ori[q_ori_ring_loc_index[[i]]])
  # }
  I_o_q_2_ori = vapply(q_ori_ring_loc_index, function(idx) mean(avg_I_2_ori[idx]), numeric(1))
  
  # Feb 9, 2026
  I_o_q_2_ori_min = min(I_o_q_2_ori)
  
  B_est_ini = max(2*abs(I_o_q_2_ori_min),10^{-20})
  A_est_ini = 2*(I_o_q_2_ori - B_est_ini/2)
  
  # Mar 14, 2026
  q_reach_thr = which(cumsum(A_est_ini) / sum(A_est_ini) >= threshold)[1]
  
  tail_idx_select = which(A_est_ini[(q_reach_thr+1):len_q] >= 0.001) 
  
  if (length(tail_idx_select) == 0) {
    num_q_max = q_reach_thr
  } else {
    num_q_max = q_reach_thr + max(tail_idx_select)
  }
  
  # tail_start = max(ceiling(0.9 * len_q), q_reach_thr + 1)
  # 
  # if (tail_start > len_q) {
  #   num_q_max = q_reach_thr
  # } else {
  #   tail_idx_select = which(A_est_ini[tail_start:len_q] >= 0.001) 
  #   
  #   if (length(tail_idx_select) == 0) {
  #     num_q_max = q_reach_thr
  #   } else {
  #     num_q_max = tail_start - 1 + max(tail_idx_select)
  #   }
  # }
  
  
  if(num_q_max/len_q<=beta){
    num_q_max=ceiling(beta*len_q)
  }
  
  ##get unique index, this maybe improved
  q_ori_ring_loc_unique_index=as.list(1:len_q)
  for(i in 1:len_q){
    unique_val=unique(avg_I_2_ori[q_ori_ring_loc_index[[i]]])
    unique_val=unique_val[1:(length(q_ori_ring_loc_index[[i]])/2)]
    q_ori_ring_loc_unique_index[[i]] = match(unique_val, avg_I_2_ori)
    # index_selected=NULL
    # for(j in 1:length(unique_val)){
    #   index_selected=c(index_selected,which(avg_I_2_ori==unique_val[j])[1])
    # }
    # q_ori_ring_loc_unique_index[[i]]=index_selected ###q_ori_ring_loc_index[[i]][index_selected]
  }
  
  total_q_ori_ring_loc_unique_index=NULL
  for(i in 1:len_q){
    total_q_ori_ring_loc_unique_index=c(total_q_ori_ring_loc_unique_index, q_ori_ring_loc_unique_index[[i]])
  }
  
  ans_list=list()
  ans_list$A_est_ini=A_est_ini
  ans_list$B_est_ini=B_est_ini
  ans_list$num_q_max=num_q_max
  ans_list$I_o_q_2_ori=I_o_q_2_ori
  ans_list$avg_I_2_ori=avg_I_2_ori
  ans_list$q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index
  ans_list$total_q_ori_ring_loc_unique_index=total_q_ori_ring_loc_unique_index
  
  #for anisotropic
  if(anisotropic){
    #image2D(t(q_ori_ring_loc[1:len_q,1:len_q]))
    q1_unique_index=as.list(1:len_q)
    q2_unique_index=as.list(1:len_q)
    
    #sub_q_ori_ring_loc=q_ori_ring_loc[1:len_q,1:len_q]
    for(i in 1:len_q){
      #index_here=which(sub_q_ori_ring_loc==i)
      #col_here=floor(index_here/len_q)+1
      #row_here=index_here%%len_q
      index_here=q_ori_ring_loc_unique_index[[i]]
      total_num_unique_index_here=length( q_ori_ring_loc_unique_index[[i]])
      q1_unique_index[[i]]=q2_unique_index[[i]]=rep(NA,total_num_unique_index_here)
      for(j in 1:(total_num_unique_index_here)){
        q2_unique_index[[i]][j]=floor((index_here[j]-1)/sz) ##could contain zero
        
        left_here=(index_here[j]-1)%%sz
        if(left_here<=len_q){
          q1_unique_index[[i]][j]=(index_here[j]-1)%%sz ##could contain zero
        }else{
          q1_unique_index[[i]][j]=sz-(index_here[j]-1)%%sz-1 ##could contain zero
        }
        # if(j<=length(row_here)){
        #   q1_unique_index[[i]][j]=row_here[j]
        #   q2_unique_index[[i]][j]=col_here[j]
        # }else{
        #   j_here=2*length(row_here)-j ###length(row_here)-(j-)
        #   q1_unique_index[[i]][j]=row_here[j_here]
        #   q2_unique_index[[i]][j]=col_here[j_here]
        # }
      }
    }
    ans_list$q1_unique_index=q1_unique_index
    ans_list$q2_unique_index=q2_unique_index
    
  }
  
  
  
  return(ans_list)
}

Get_MSD_nonparametric_mf<-function(theta,d_input,input_training,model_name){
  if(model_name=='direct_nonparametric'){
    ##range.par can be adjust
    #100*d_input[2]
    range_par=(max(input_training)-min(input_training))
    
    ##trend
    #X=cbind(rep(1,length(input_training)),as.numeric(input_training))
    #X_testing=cbind(rep(1,length(d_input[-1])),as.numeric(log(d_input[-1])))
    #
    #m_direct=rgasp(design=input_training,response=log(theta),
    #               range.par = range_par,trend=as.matrix(X)) ##input_training is log
    #m_direct_log_MSD=predict(m_direct,testing_input=as.matrix(log(d_input[-1])),testing_trend=X_testing)$mean
    
    ##no trend
    m_direct=RobustGaSP::rgasp(design=input_training,response=log(theta),range.par = range_par) ##input_training is log 
    
    ##small nugget
    #m_direct=rgasp(design=input_training,response=log(theta),range.par = range_par,nugget=10^{-4}) ##input_training is log
    
    m_direct_log_MSD=RobustGaSP::predict(m_direct,as.matrix(log(d_input[-1])))$mean
    
    MSD=c(0,exp(m_direct_log_MSD))
    
  }
  # else if(model_name=='direct_additive_nonparametric'){
  #   ##range.par can be adjust
  #   #100*d_input[2]
  #   range_par=(max(input_training)-min(input_training))
  #   output=log(cumsum(theta))
  #   m_direct=rgasp(design=input_training,response=output,range.par = range_par) ##input_training is log
  #   m_direct_log_MSD=predict(m_direct,as.matrix(log(d_input[-1])))$mean
  #   MSD=c(0,exp(m_direct_log_MSD))
  # }
  else if(model_name=='varying_power_nonparametric'){
    #if(input_training[1]==0){
    ###this is just for time to be 0 to 1, we can always input 0 to 1 and transform later
    range_par=(max(input_training)-min(input_training)) ##input training is in the log space
    
    alpha_training=2*theta[-1]/(1+theta[-1])
    m_varying_power=RobustGaSP::rgasp(design=input_training[-1],response=(alpha_training),range.par = range_par) ##input_training is log
    m_varying_power_pred_alpha=RobustGaSP::predict(m_varying_power,as.matrix(log(d_input[-1])))$mean
    log_msd=log(theta[1])+(m_varying_power_pred_alpha)*log(d_input[-1])
    
    MSD=c(0,exp(log_msd))
    #}else{ ##need to make sure no 0 is in training
    
    # }
    
  }
  else if(model_name=='varying_power_nonparametric_msd'){
    range_par=(max(input_training)-min(input_training)) ##input training is in the log space
    log_MSD_training=rep(NA,length(input_training))
    log_MSD_training[1]=log(theta[1])
    
    log_MSD_training[2:length(input_training)]= log_MSD_training[1]+2*theta[-1]/(1+theta[-1])*(input_training[-1])
    
    m_direct=RobustGaSP::rgasp(design=input_training,response=(log_MSD_training),range.par = range_par) ##input_training is log
    m_direct_log_MSD=RobustGaSP::predict(m_direct,as.matrix(log(d_input[-1])))$mean
    
    #X=cbind(rep(1,length(input_training)),input_training)
    #X_testing=cbind(rep(1,length(d_input[-1])),log(d_input[-1]))
    
    
    #m_direct=rgasp(design=input_training,response=(log_MSD_training),
    #               range.par = range_par*0.1,trend=X) ##input_training is log
    #m_direct_log_MSD=predict(m_direct,as.matrix(log(d_input[-1])),testing_trend=X_testing)$mean
    
    MSD=c(0,exp(m_direct_log_MSD))
    
  }
  return(MSD)
}


# get_circ_firstrow = function(p, first_row){
#   second_half = rev(first_row[2:p])
#   circ_firstrow = c(first_row, second_half)
#   return(circ_firstrow)
# }

###changed July 15, 2025
log_lik_param_nonparametric_ori_mf<-function(param,I_q_cur,B_cur,num_q_cur,I_o_q_2_ori,
                                          d_input,q_ori_ring_loc_unique_index,sz,len_t,q,model_name,input_training,
                                          q_index_selected=1:num_q_cur,subsample_t=1,whitening=F
){ #recompute_num_q_cur=F
  #print(param)
  #theta=exp(param)
  p=length(param)-1
  #theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
  
  if(whitening){
    range_par=(max(input_training)-min(input_training))
    ##no trend
    m_direct=RobustGaSP::rgasp(design=input_training,response=param[-(p+1)],range.par = range_par) ##input_training is log
    log_theta=m_direct@L%*%param[-(p+1)]
    theta=exp(log_theta)
  }else{
    theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
    
  }
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }
  
  A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
  
  
  
  ##the model is defined by MSD
  #MSD=Get_MSD_nonparametric_mf(theta,d_input,input_training,model_name)
  
  subsample_d_index=seq(1,length(d_input),subsample_t)
  d_input_subsample=d_input[subsample_d_index]
  
  MSD=Get_MSD_nonparametric_mf(theta,d_input_subsample,input_training,model_name)
  
  
  log_lik_sum=0
  len_t_subsample=length(MSD)
  
  NTz <- SuperGauss::NormalToeplitz$new(len_t_subsample) ##maybe create it outside?
  
  eta=B_cur/4 ##nugget
  
  num_q_select=length(q_index_selected)
  
  i_log_lik_record=rep(NA,num_q_select)
  
  for(i_q_selected in q_index_selected){
    #print(q_selected)
    
    #output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
    #output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
    output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],subsample_d_index])/(sz)
    output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],subsample_d_index])/(sz)
    
    
    q_selected=q[i_q_selected]
    #beta_q = (D*q[i_q_selected]^2)
    sigma_2=A_cur[i_q_selected]/4
    
    acf = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
    acf[1] = acf[1]+eta
    acf=as.numeric(acf)
    
    i_log_lik_record[i_q_selected] = sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
    
    # if (is.nan(i_log_lik)) {
    #   ## try new estimator: get the first row of the new estimator
    #   circ_mat = generate_circ_mat(length(acf), acf)
    #   lambda = fft(circ_mat[1,]) ## eigen values
    #   lambda_truncated = ifelse(Re(lambda) <= 0, 10^{-4}, lambda)
    #   firstrow_reconstructed = Re(fft(lambda_truncated, inverse = TRUE) / (2*length(acf) - 1)) # the first row of reconstructed circulant matrix
    #   i_log_lik = sum(NTz$logdens(z = t(output_re), acf = firstrow_reconstructed[1:length(acf)])) + sum(NTz$logdens(z = t(output_im), acf = firstrow_reconstructed[1:length(acf)]))
    # }
    
    #log_lik_sum=log_lik_sum+sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
    log_lik_sum = log_lik_sum +  i_log_lik_record[i_q_selected]
    
  }
  
  #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index))*log(2*pi) ##add 2pi
  #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index)*len_t)*log(2*pi) ##add 2pi? we may not need to add these constant, we should remove them later
  if(is.nan(log_lik_sum)){
    log_lik_sum=-10^15
    #log_lik_sum=-10^15-10^15*runif(1)
    
    
  }
  
  #print(c((param),log_lik_sum))
  #print(c(log_lik_sum))
  return(log_lik_sum)
}



Get_loglik_nonparametric_mf <- function(I_q_cur, B_cur, num_q_cur, I_o_q_2_ori,
                                     d_input, q_ori_ring_loc_unique_index, sz, len_t, q,
                                     model_name, input_training, q_index_selected,
                                     subsample_t = 1, substack = FALSE, whitening = FALSE
) {
  if (isTRUE(substack)) {
    subsample_d_index = seq(1, length(d_input) - 1, subsample_t)
    d_input_subsample = d_input[c(1, subsample_d_index + 1)]
    time_index = c(1, subsample_d_index + 1)
  } else {
    subsample_d_index = seq(1, length(d_input), subsample_t)
    d_input_subsample = d_input[subsample_d_index]
    time_index = subsample_d_index
  }
  
  len_t_subsample = length(d_input_subsample)
  NTz = SuperGauss::NormalToeplitz$new(len_t_subsample)
  
  re_list = vector("list", length(q_index_selected))
  im_list = vector("list", length(q_index_selected))
  q_vals = q[q_index_selected]
  
  for (ii in seq_along(q_index_selected)) {
    iq = q_index_selected[ii]
    idx = q_ori_ring_loc_unique_index[[iq]]
    re_list[[ii]] = Re(I_q_cur[idx, time_index]) / sz
    im_list[[ii]] = Im(I_q_cur[idx, time_index]) / sz
  }
  
  function(param) {
    p = length(param) - 1
    
    if (isTRUE(whitening)) {
      range_par = max(input_training) - min(input_training)
      m_direct = RobustGaSP::rgasp(
        design = input_training,
        response = param[-(p + 1)],
        range.par = range_par
      )
      log_theta = m_direct@L %*% param[-(p + 1)]
      theta = exp(log_theta)
    } else {
      theta = exp(param[-(p + 1)])
    }
    
    B_here = B_cur
    if (is.na(B_here)) {
      sigma_2_0_hat = exp(param[p + 1])
      B_here = 2 * sigma_2_0_hat
    }
    
    A_cur = abs(2 * (I_o_q_2_ori - B_here / 2))
    MSD = Get_MSD_nonparametric_mf(theta, d_input_subsample, input_training, model_name)
    
    log_lik_sum = 0
    
    for (ii in seq_along(q_index_selected)) {
      iq = q_index_selected[ii]
      sigma_2 = A_cur[iq] / 4
      acf = sigma_2 * exp(-q_vals[ii]^2 * MSD / 4)
      acf[1] = acf[1] + B_here / 4
      acf = as.numeric(acf)
      
      log_lik_sum = log_lik_sum +
        sum(NTz$logdens(z = t(re_list[[ii]]), acf = acf)) +
        sum(NTz$logdens(z = t(im_list[[ii]]), acf = acf))
    }
    
    if (is.nan(log_lik_sum)) {
      log_lik_sum = -1e15
    }
    
    log_lik_sum
  }
}



get_initial_param_nonparametric_mf <- function(model_name,d_input,num_estim=6,sigma_0_2_ini=NA,design_optimization="log_equal_space"){
  if(design_optimization=='equal_space'){
    input_index=rep(NA,num_estim)
    input_index[1]=2 ##perhaps start from the beginning, this is like half or full
    d_int=floor(length(d_input)/(num_estim-1))
    input_index[2:num_estim]=floor(1:(num_estim-1))*d_int
  }else if(design_optimization=='log_equal_space'){ ##log equal space 
    input_index = 1 + unique(ceiling(exp(seq(from=log(1),to=log(length(d_input)-1),
                                             by=log(length(d_input)-1)/(num_estim-1))))) 
    
    input_index[length(input_index)] = length(d_input)
    
  }else if(design_optimization=='log_equal_space_middle'){ ##middle 
    skip=8
    d_int=(log(d_input[length(d_input)])-log(d_input[2]))/(num_estim-1+skip)
    
    input_index=rep(NA,num_estim)
    
    input_index[1]=2 ##perhaps start from the beginning, first one
    input_seq_fit=log(d_input[2])+d_int*(skip:(num_estim+skip-1)) ##the subsetting is like skip some in the log space? this is almost like equally space then, perhaps try equally space
    for(i in 2:(length(input_index))){
      input_index[i]=which(abs(log(d_input)-input_seq_fit[i-1])==min(abs(log(d_input)-input_seq_fit[i-1])))
    }
  }
  input_training=log(d_input[input_index])
  
  
  if(model_name=='direct_nonparametric'){
    param_initial=log(c((exp(input_training)*0.1),sigma_0_2_ini)) ##
  }else if(model_name=='varying_power_nonparametric'){
    param_initial1=log(c(0.1,rep(.1,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
    param_initial2=log(c(0.1,rep(1,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
    
    param_initial=cbind(param_initial1,param_initial2)
  }else if(model_name=='varying_power_msd_nonparametric'){
    param_initial=log(c(0.1,rep(1,length(input_training)-1),sigma_0_2_ini)) ##ad hoc
    
  }
  return.list=list()
  return.list$input_training=input_training
  return.list$param_initial=param_initial
  return.list$input_index=input_index
  return(return.list)
}



param_uncertainty_nonparametric_mf <- function(param_est,I_q_cur,B_cur=NA,num_q_cur,I_o_q_2_ori,
                                          q_ori_ring_loc_unique_index,
                                          sz,len_t,d_input,q,
                                          q_index_selected, 
                                          model_name,input_training,estimation_method='asymptotics',M,num_iteration_max,lower_bound=NA){
  
  
  
  q_lower=q-min(q)
  
  loglik_q_lower = Get_loglik_nonparametric_mf(I_q_cur = I_q_cur, B_cur = NA, num_q_cur = num_q_cur,
                                            I_o_q_2_ori = I_o_q_2_ori, q_ori_ring_loc_unique_index = q_ori_ring_loc_unique_index,
                                            sz = sz, len_t = len_t, d_input = d_input, q = q_lower,
                                            model_name = "direct_nonparametric", input_training = input_training,
                                            q_index_selected = q_index_selected
  )
  
  m_param_lower = optim(param_est, loglik_q_lower, control = list(fnscale = -1, maxit = 100), method = "L-BFGS-B")
  
  q_upper=q+min(q)
  
  loglik_q_upper = Get_loglik_nonparametric_mf(I_q_cur = I_q_cur, B_cur = NA, num_q_cur = num_q_cur,
                                            I_o_q_2_ori = I_o_q_2_ori, q_ori_ring_loc_unique_index = q_ori_ring_loc_unique_index,
                                            sz = sz, len_t = len_t, d_input = d_input, q = q_upper,
                                            model_name = "direct_nonparametric", input_training = input_training,
                                            q_index_selected = q_index_selected
  )
  
  m_param_upper = optim(param_est, loglik_q_lower, control = list(fnscale = -1, maxit = 100), method = "L-BFGS-B")
  
  
  ##here this is log theta 
  p=length(param_est)-1
  param_range=matrix(NA,2,p+1)
  for(i in 1:(p+1) ){
    param_range[1,i]=min(m_param_lower$par[i],m_param_upper$par[i])
    param_range[2,i]=max(m_param_lower$par[i],m_param_upper$par[i])
    
  }
  half_length_param_range_fft=(param_range[2,]-param_range[1,])/2
  
  ##this asymptotics is on log theta as well
  if(estimation_method=='asymptotics'){
    
    p=length(param_est)-1
    theta=exp(param_est[-(p+1)]) ##first p parameters are parameters in ISF
    if(is.na(B_cur)){ ##this fix the dimension
      sigma_2_0_hat=exp(param_est[p+1]) ##noise
      B_cur=2*sigma_2_0_hat
    }
    #A_cur = 2*(I_o_q_2_ori - B_cur/2)
    
    A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
    eta=B_cur/4 ##nugget
    
    #why rescaled? only MSD_grad matters, MSD and grad_trans do not matter
    #MSD=Get_MSD_nonparametric(theta,d_input=d_input_rescale,model_name=model_name,input_training=input_training_rescale)
    #MSD_grad_save=Get_MSD_grad_nonparameteric(theta,d_input_rescale,model_name,input_training=input_training_rescale)
    #grad_trans_save=Get_grad_trans_nonparametric(theta,d_input_rescale,model_name)
    grad_trans=rep(1,p) ##derivative on the param not theta
    MSD=Get_MSD_nonparametric_mf(theta,d_input=d_input,model_name=model_name,input_training=input_training)
    ##numerical gradient
    delta=0.001 ## delta for numerical difference
    MSD_grad=matrix(NA,length(MSD),p)
    for(i_p in 1:p){
      param_change=param_est
      param_change[i_p]=  param_change[i_p]+delta
      theta_change=exp(param_change[-(p+1)]) ##first p parameters are parameters in ISF
      # if((B_cur==2*sigma_2_0_hat)){ ##this fix the dimension
      #   sigma_2_0_change=exp(param_change[p+1]) ##noise
      #   B_change=2*sigma_2_0_change
      # }
      MSD_change=Get_MSD_nonparametric_mf(theta_change,d_input=d_input,model_name=model_name,input_training=input_training)
      MSD_grad[,i_p]=(MSD_change-MSD)/delta
      
    }
    #MSD=Get_MSD_nonparametric(theta,d_input=d_input,model_name=model_name,input_training=input_training)
    #MSD_grad=Get_MSD_grad_nonparameteric(theta,d_input,model_name,input_training=input_training)
    #grad_trans=Get_grad_trans_nonparametric(theta,d_input,model_name)
    
    #Hessian_inv_list=as.list(1:num_q_cur)
    #Hessian_inv_sum=0
    Hessian_list=as.list(1:num_q_cur)
    Hessian_sum=0
    
    for(i_q_selected in 1: num_q_cur){
      #print(q_selected)
      
      #output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
      #output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
      
      q_selected=q[i_q_selected]
      
      sigma_2=A_cur[i_q_selected]/4
      
      acf0 = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
      acf=acf0
      acf[1] = acf[1]+eta ##for grad this is probably no adding
      acf=as.numeric(acf)
      
      Tz <- SuperGauss::Toeplitz$new(len_t,acf=acf) ##maybe create it outside?
      Hessian=matrix(NA,p+1,p+1) ##last one is
      acf_grad=matrix(NA,len_t,p+1)
      for(i_p in 1:p){
        acf_grad[,i_p]=-acf0*q_selected^2/4* MSD_grad[,i_p]*grad_trans[i_p]
      }
      #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
      acf_grad[,p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2)
      #acf_grad[,p+1]=(-0.5)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2)
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
      
      
      #Hessian_inv_list[[i_q_selected]]=solve(Hessian)
      #Hessian_inv_sum=Hessian_inv_sum+ Hessian_inv_list[[i_q_selected]]*length(q_ori_ring_loc_unique_index[[i_q_selected]])
      ##a litte more conservation is to say they are perfectly correlated in a ring
      #Hessian_inv_sum=Hessian_inv_sum+ Hessian_inv_list[[i_q_selected]]
      
      
      #cor(Re(I_q_cur[q_ori_ring_loc_unique_index[[1]][2],]),Re(I_q_cur[q_ori_ring_loc_unique_index[[1]][1],]))
    }
    
    Hessian_sum=Hessian_sum*M/sum(lengths(q_ori_ring_loc_unique_index[1:num_q_cur]))
    sd_theta_B=sqrt(diag(solve(Hessian_sum)))
    #sd_theta_B=sqrt(diag(Hessian_inv_sum/sum(lengths(q_ori_ring_loc_unique_index))^2 ))
    
    param_range_min_max = matrix(nrow = 2, ncol = ncol(param_range))
    param_range_w_asymp = matrix(nrow = 2, ncol = ncol(param_range))
    for(i_p in 1:length(param_est)){
      param_range_min_max[,i_p]=c(min(param_est[i_p], param_range[1,i_p]),
                                  max(param_est[i_p], param_range[2,i_p]))
    }
    
    param_range_w_asymp[1,]=param_range_min_max[1,]-sd_theta_B*qnorm(0.975) 
    param_range_w_asymp[2,]=param_range_min_max[2,]+sd_theta_B*qnorm(0.975)
    
    theta_lower=exp(param_range_w_asymp[1,-ncol(param_range)])
    theta_upper=exp(param_range_w_asymp[2,-ncol(param_range)])
    
    MSD_lower=Get_MSD_nonparametric_mf(theta=theta_lower,d_input=d_input,input_training=input_training,model_name=model_name)
    MSD_upper=Get_MSD_nonparametric_mf(theta=theta_upper,d_input=d_input,input_training=input_training,model_name=model_name)
    
    for (i in 1:length(MSD)) {
      MSD_lower[i] = min(MSD_lower[i], MSD[i])
      MSD_upper[i] = max(MSD_upper[i], MSD[i])
    }
  }
  
  result = list()
  result$param_range = param_range_w_asymp
  result$MSD_lower = MSD_lower
  result$MSD_upper = MSD_upper
  return(result)
}



matern_5_2 <- function(xi, xj, sigma, l) {
  d  = abs(outer(xi, xj, "-"))
  r  = d / l
  K  = sigma^2 * (1 + sqrt(5)*r + (5/3)*r^2) * exp(-sqrt(5)*r)
  return(K)
}

## covariance matrix for log MSD
build_cov_logMSD <- function(x, sd_log_MSD, l) {
  Corr_mat = matern_5_2(x, x, sigma=1, l=l)
  K = Corr_mat * outer(sd_log_MSD, sd_log_MSD)
  K = K + 1e-8 * diag(length(x))
  return(K)
}


# sampling log MSD
sample_logMSD_matern <- function(y, x, sd_log_MSD, l, n_samples) {
  n = length(x)
  
  # Covariance matrix
  K = build_cov_logMSD(x, sd_log_MSD=sd_log_MSD, l=l)
  L = t(chol(K))
  
  samples = matrix(0, nrow=n_samples, ncol=n)
  for (i in seq_len(n_samples)) {
    z = rnorm(n)
    samples[i, ] = y + L %*% z
  }
  
  return(samples)
}


calculate_Gp_Gpp_GSER_from_MSD <- function(MSD_est, d_input, MSD_unit = 1,
                                          kB = 1.380649*10^{-23},
                                          Temp = 1, a = 1) {
  MSD_est = MSD_est * MSD_unit
  idx = which(MSD_est > 0)
  n = length(idx)
  
  MSD_ln_dev_est = numeric(n)
  ln_MSD_est = log(MSD_est[idx])
  ln_time = log(d_input[idx])
  
  if (n >= 3) {
    MSD_ln_dev_est[2:(n - 1)] = (ln_MSD_est[3:n] - ln_MSD_est[1:(n - 2)])/(ln_time[3:n] - ln_time[1:(n - 2)])
  }
    
  MSD_ln_dev_est[1] = (ln_MSD_est[2] - ln_MSD_est[1])/(ln_time[2] - ln_time[1])
  MSD_ln_dev_est[n] = (ln_MSD_est[n] - ln_MSD_est[n - 1])/(ln_time[n] - ln_time[n - 1])
  
  
  Gcomplex_est = 2 * kB * Temp/(3 * pi * a * MSD_est[idx] * gamma(1 + MSD_ln_dev_est))
  
  Gprime_est = Gcomplex_est * cos(pi * MSD_ln_dev_est / 2)
  Gdoubleprime_est = Gcomplex_est * sin(pi * MSD_ln_dev_est / 2)
  omega = 2 * pi / d_input[idx]
  
  return_list = list()
  return_list$omega = omega
  return_list$Gprime_est = Gprime_est
  return_list$Gdoubleprime_est = Gdoubleprime_est
  
  return_list
}


sample_Gp_Gpp_GSER <- function(y, x, d_input, sd_log_MSD, l, n_samples, MSD_unit = 1, kB = 1.380649*10^{-23},
                              Temp = 1, a = 1) {
  logMSD_samples = sample_logMSD_matern(y, x, sd_log_MSD, l, n_samples)
  MSD_samples = exp(logMSD_samples)
  
  Gp_list = matrix(NA, nrow = n_samples, ncol = ncol(logMSD_samples))
  Gpp_list = matrix(NA, nrow = n_samples, ncol = ncol(logMSD_samples))
  omega = NULL
  
  for (m in seq_len(n_samples)) {
    res = calculate_Gp_Gpp_GSER_from_MSD(MSD_est = MSD_samples[m, ], d_input = d_input, MSD_unit = MSD_unit,
                                         kB = kB, Temp = Temp, a = a)
    
    omega = res$omega
    Gp_list[m, seq_along(res$Gprime_est)] = res$Gprime_est
    Gpp_list[m, seq_along(res$Gdoubleprime_est)] = res$Gdoubleprime_est
  }
  
  return_list = list()
  return_list$omega = omega
  return_list$Gp = Gp_list
  return_list$Gpp = Gpp_list
  
  return_list
}


#' Estimate Viscoelastic Moduli by GSER
#'
#' @description
#' Estimate the storage and loss moduli from mean squared displacement (MSD)
#' using the generalized Stokes-Einstein relation. Because the transformation
#' from MSD to the storage and loss moduli is nonlinear, log MSD curves are samples
#' and used to calculate storage and loss moduli. The pointwise medians of the sampled 
#' moduli are used as the estimates.
#'
#' @param MSD_est A numeric vector of estimated MSD values. This argument is used
#'   for compatibility and for computing theoretical moduli when a simulation
#'   object is provided. The estimated moduli are obtained from
#'   \code{MSD_lower} and \code{MSD_upper}.
#' @param d_input A numeric vector of lag times.
#' @param MSD_lower A numeric vector of lower bound of 95\% confidence interval of MSD.
#' @param MSD_upper A numeric vector of upper bound of 95\% confidence interval of MSD.
#' @param sim_object Optional simulation object. If provided, the function also
#'   returns the theoretical storage and loss moduli for the simulated data.
#' @param MSD_unit Scaling factor applied to MSD values used to convert the
#'   estimated MSD to the desired physical unit before modulus estimation.
#' @param kB Boltzmann constant.
#' @param Temp Absolute temperature.
#' @param a Particle radius.
#'
#' @return
#' A list containing estimated frequency, storage modulus, and loss modulus.
#' If \code{sim_object} is provided, the theoretical moduli are also returned.
#'
#' @export


Get_Gp_Gpp_GSER <- function(MSD_est, d_input, MSD_lower, MSD_upper,
                           sim_object = NA, MSD_unit = 1,
                           kB = 1.380649*10^{-23}, Temp = 1, a = 1) {
  
  if (missing(MSD_lower) || missing(MSD_upper)) {
    stop("MSD_lower and MSD_upper must be provided for GSER modulus estimation.")
  }
  
  set.seed(1)
  n_samples = 1000
  
  idx = which(d_input > 0 & MSD_lower > 0 & MSD_upper > 0)
  x = log(d_input[idx])
  y = (log(MSD_upper[idx]) + log(MSD_lower[idx])) / 2
  sd_log_MSD = (log(MSD_upper[idx]) - log(MSD_lower[idx])) / (2 * 1.96)
  l_range = 0.1 * diff(range(x))
  
  samples_res = sample_Gp_Gpp_GSER(y = y, x = x, d_input = d_input[idx], sd_log_MSD = sd_log_MSD,
                                   l = l_range, n_samples = n_samples, MSD_unit = MSD_unit,
                                   kB = kB, Temp = Temp, a = a)
  
  Gp_samples = samples_res$Gp
  Gpp_samples = samples_res$Gpp
  
  Gp_samples[Gp_samples <= 0] = NA
  Gpp_samples[Gpp_samples <= 0] = NA
  
  return_list = list()
  return_list$omega = samples_res$omega
  return_list$Gprime_est = apply(Gp_samples, 2, stats::median, na.rm = TRUE)
  return_list$Gdoubleprime_est = apply(Gpp_samples, 2, stats::median, na.rm = TRUE)
  
  if (methods::is(sim_object, "simulation")) {
    n = length(idx)
    simul_model_name = sim_object@model_name
    
    if (simul_model_name == "BM") {
      sigma_bm = sim_object@param[1]
      beta = 2 * sigma_bm^2
      
      MSD_truth = beta * d_input[idx]
      MSD_ln_dev_truth = 1
    } else if (simul_model_name == "FBM") {
      sigma_fbm = sim_object@param[1]
      H = sim_object@param[2]
      
      beta = 2 * sigma_fbm^2
      alpha = 2 * H
      
      MSD_truth = beta * d_input[idx]^alpha
      MSD_ln_dev_truth = alpha
    } else if (simul_model_name == "OU") {
      rho = sim_object@param[1]
      sigma_ou = sim_object@param[2]
      amplitude = 4 * sigma_ou^2
      
      MSD_truth = amplitude * (1 - rho^d_input[idx])
      MSD_ln_dev_truth = d_input[idx] * (-log(rho)) * rho^d_input[idx] /
        (1 - rho^d_input[idx])
    } else if (simul_model_name == "OU+FBM") {
      rho = sim_object@param[1]
      sigma_ou = sim_object@param[2]
      sigma_fbm = sim_object@param[3]
      H = sim_object@param[4]
      
      amplitude = 4 * sigma_ou^2
      beta = 2 * sigma_fbm^2
      alpha = 2 * H
      
      MSD_truth = beta * d_input[idx]^alpha + amplitude * (1 - rho^d_input[idx])
      dMSD_d_input = beta * alpha * d_input[idx]^(alpha - 1) + amplitude * (-log(rho)) * rho^d_input[idx]
      MSD_ln_dev_truth = (d_input[idx] / MSD_truth) * dMSD_d_input
    }
    
    Gcomplex_truth = 2 * kB * Temp / (3 * pi * a * MSD_truth * gamma(1 + MSD_ln_dev_truth))
    Gprime_truth = Gcomplex_truth * cos(pi * MSD_ln_dev_truth / 2)
    Gdoubleprime_truth = Gcomplex_truth * sin(pi * MSD_ln_dev_truth / 2)
    
    if (simul_model_name == "BM") {
      Gprime_truth = rep(0, n)
      return_list$Gprime_est = rep(0, length(return_list$omega))
    }
    
    return_list$Gprime_truth = Gprime_truth
    return_list$Gdoubleprime_truth = Gdoubleprime_truth
  }
  
  return_list
}





Get_initial_param_DDM <- function(params = list()) {
  
  DDM_result = list()
  
  if (is.null(params$pxsz)) DDM_result$pxsz = 1 else DDM_result$pxsz = params$pxsz
  if (!is.null(params$nframes)) DDM_result$nframes = params$nframes else stop("Please specify number of time frames.")
  if (!is.null(params$nx)) DDM_result$nx = params$nx else stop("Please specify size.")
  if (!is.null(params$ny)) DDM_result$ny = params$ny else stop("Please specify size.")
  if (!is.null(params$mindt)) DDM_result$mindt = params$mindt else stop("Please specify mininum lag time.")
  if (!is.null(params$q)) DDM_result$q = params$q else stop("Please specify q.")
  if (!is.null(params$len_q)) DDM_result$len_q = params$len_q else stop("Please specify len_q.")
  if (!is.null(params$I_q_matrix)) DDM_result$I_q_matrix = params$I_q_matrix else stop("Please specify I_q_matrix.")
  if (!is.null(params$q_ori_ring_loc_index)) DDM_result$q_ori_ring_loc_index = params$q_ori_ring_loc_index else stop("Please specify q_ori_ring_loc_index.")
  if (!is.null(params$index_dt_selected)) DDM_result$index_dt_selected = params$index_dt_selected else DDM_result$index_dt_selected = 1:2 
  if (!is.null(params$I_o_q_2)) DDM_result$I_o_q_2 = params$I_o_q_2 
  
  
  stopifnot(is.numeric(DDM_result$pxsz), is.numeric(DDM_result$mindt), 
            is.numeric(DDM_result$nframes), is.numeric(DDM_result$nx), is.numeric(DDM_result$ny))
  
  DDM_result$dt = DDM_result$mindt * (1:(DDM_result$nframes - 1))
  DDM_result$ndt = (DDM_result$nframes - 1):1
  
  if (min(DDM_result$nx, DDM_result$ny) %% 2 == 0) {
    DDM_result$sz = min(DDM_result$nx, DDM_result$ny) - 1
  } else {
    DDM_result$sz = min(DDM_result$nx, DDM_result$ny)
  }
  
  
  DDM_result$Dqt = matrix(nrow = DDM_result$len_q, ncol = length(DDM_result$index_dt_selected))

  for (idx in seq_along(DDM_result$index_dt_selected)) {
    k = DDM_result$index_dt_selected[idx]
    
    avg = numeric(nrow(DDM_result$I_q_matrix))
    for (j in seq_len(DDM_result$ndt[k])) {
      diff = DDM_result$I_q_matrix[, j + k] - DDM_result$I_q_matrix[, j]
      avg = avg + abs(diff)^2
    }
    avg = avg / (DDM_result$ndt[k] * DDM_result$sz^2)
    
    DDM_result$Dqt[, idx] = vapply(DDM_result$q_ori_ring_loc_index, function(ii) mean(avg[ii]), numeric(1))
  }
  
  ddm_ori = DDM_result$Dqt
  q_ori = DDM_result$q
  I_o_2_ori = matrix(DDM_result$I_o_q_2, nrow=1)
  
  # A_hat, B_hat 
  DDM_result$B_hat = min(mean(ddm_ori[nrow(ddm_ori),]), min(ddm_ori[,1]))
  DDM_result$sigma_2_0_hat = DDM_result$B_hat/2
  DDM_result$A_hat = 2 * (I_o_2_ori - DDM_result$sigma_2_0_hat)
  
  n_q = length(q_ori)
  n_q_min = 4
  n_q_max = min(n_q, 83)
  Dqt_select_dt = ddm_ori[n_q_min:n_q_max, DDM_result$index_dt_selected, drop=FALSE]
  A_hat = as.numeric(DDM_result$A_hat)[n_q_min:n_q_max]
  q = q_ori[n_q_min:n_q_max]
  
  # f(q,t)
  fqt = matrix(nrow = nrow(Dqt_select_dt), ncol = ncol(Dqt_select_dt))
  for (i in 1:nrow(Dqt_select_dt)) {
    fqt[i, ] = 1 - (Dqt_select_dt[i, ] - DDM_result$B_hat) / A_hat[i]
  }
  fqt[fqt <= 0] = NA
  
  # MSD(q,t)
  Qmat = matrix(q, nrow=length(q), ncol=length(DDM_result$index_dt_selected))
  MSD_q_dt = 4 * log(1/fqt) / (Qmat^2)
  
  MSD = rep(NA, length(DDM_result$index_dt_selected))
  for (ti in DDM_result$index_dt_selected) {
    MSD[ti] = stats::median(MSD_q_dt[, ti], na.rm=TRUE)
  }
  
  DDM_result$MSD = MSD
  return(DDM_result)
}



