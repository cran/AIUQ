## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(AIUQ)

## ----fig.height = 3.5, fig.width = 3.5----------------------------------------
set.seed(1)

sim_bm = simulation()
show(sim_bm)

## Plot simulated particle trajectory
plot_traj(sim_bm)

## ----fig.height = 3.5, fig.width = 7------------------------------------------
par(mfrow=c(1,2))
## Plot intensity profile for different frames 
plot_intensity(sim_bm@intensity, sz=sim_bm@sz) #first frame 
plot_intensity(sim_bm@intensity, sz=sim_bm@sz,frame=10, color=T) #10th frame, color image

## ----fig.height = 3.5, fig.width = 7------------------------------------------
## AIUQ method: use BM as fitted model 
sam = SAM(sim_object=sim_bm, model_name='BM')
show(sam)

par(mfrow=c(1,2))
## Plot true MSD and estimated MSD
plot_MSD(object=sam, msd_truth=sam@msd_truth) #in log10 scale

## Plot intensity in reciprocal space 
plot_I_q_1 = matrix(sam@I_q[,1], sam@sz[1],sam@sz[2]) #first frame 
plot3D::image2D(abs(fftshift(plot_I_q_1)),main="intensity in reciprocal space")


## -----------------------------------------------------------------------------
sam = SAM(sim_object=sim_bm, AIUQ_thr=c(0.99,0.6)) 
#Note: Default model_name is "BM", it's ok to not specify this argument if want to fit with BM model 
show(sam)

sam = SAM(sim_object=sim_bm, index_q_AIUQ=5:50)
show(sam)

## ----fig.height = 3.5, fig.width = 3.5----------------------------------------
set.seed(1)

## Simulation
sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
show(sim_bm)

## Plot simulated particle trajectory
plot_traj(sim_bm)

## ----fig.height = 3.5, fig.width = 7------------------------------------------
## AIUQ method: fitting using BM model with uncertainty quantification 
sam = SAM(sim_object=sim_bm, uncertainty=T)
show(sam)

par(mfrow=c(1,2))
## Plot true MSD and estimated MSD with uncertainty 
plot_MSD(object=sam, msd_truth=sam@msd_truth) #in log10 scale 
plot_MSD(object=sam, msd_truth=sam@msd_truth,log10=F) #in real scale

## ----fig.height = 3.5, fig.width = 7------------------------------------------
set.seed(1)

## Simulation
sim_ou = simulation(sigma_ou=4, model_name="OU")
show(sim_ou)

par(mfrow=c(1,2))
## Plot simulated particle trajectory
plot_traj(sim_ou)

## AIUQ method: fitting using OU model
sam_ou = SAM(sim_object=sim_ou, model_name=sim_ou@model_name) 
show(sam_ou)

## Plot true MSD and estimated MSD with uncertainty
plot_MSD(object=sam_ou, msd_truth=sam_ou@msd_truth) #in log10 scale

## -----------------------------------------------------------------------------
set.seed(1)

## Simulation
sim_bm = simulation(sz=100,len_t=100,sigma_bm=0.5)
show(sim_bm)

## User defined MSD structure: function of parameters and 
#                                          vector of lag times 
msd_fn = function(param, d_input){
  MSD = param[1]*d_input+param[2]*d_input^2
}

# show MSD and MSD gradient with a simple example 
theta = c(2,1)
d_input = 0:10
model_name = "user_defined"
MSD_list = get_MSD_with_grad(theta=theta,d_input=d_input,model_name=model_name,
                             msd_fn=msd_fn)
MSD_list$msd

## AIUQ method: fitting using user_defined model 
sam = SAM(sim_object=sim_bm, model_name=model_name, msd_fn=msd_fn, num_param=2)
show(sam)

