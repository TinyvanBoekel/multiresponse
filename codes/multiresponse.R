# R script belonging to the manuscript "Multiresponse kinetics with estimation of the experimental variance-covariance matrix. Full Bayesian analysis using Stan"
library(tidyverse)
library(patchwork) # to combine ggplot results
library(cmdstanr)
library(here)  # to make references to files relative
library(GGally)
library(tidybayes)
library(kableExtra)
library(bayesplot)
library(posterior)  # library for posterior analysis
library(deSolve)    # library for solving differential equations
library(gslnls)     # library for least squares fit using levenberg-marquart algorithm
# Load Box-Draper data set:
BD_df <- read.csv(file=here("data","ABC.csv"), header=TRUE, sep=",") |> as_tibble()
# format neeeded to do regression with gslnls:
BD_obs <- pivot_longer(BD_df, cols = c("A", "B", "C"), names_to = "species", cols_vary = "slowest")
# make plot of Box-Draper data set
BD_obs %>% ggplot(aes(x=time, y=value, colour=species, shape=species))+
  geom_point(size=2)+
  labs(x="time (arbitrary values)", y="[A],[B],[C] (arbitrary values)")
# compile Stan file:
BD_model <- cmdstan_model(here("fits","BD.stan"))
# Prepare data for Stan to sample:
t_pred = seq(from =0, to = max(BD_df$time), length.out=40)
BD_data <- list(
  N_t = length(BD_df$time),
  N_r=3,
  ts = BD_df$time,
  y = BD_df[,2:4],
  y0=1.0,
  t_pred = t_pred,
  N_p=length(t_pred),
  lkj_df=2
)
# sampling from the compiled model:
BD_model_regr <- BD_model$sample(
  data = BD_data,
  seed = 42L,
  refresh = 0,                # silence progress notification
  parallel_chains=4 
)
#Store regression model:
BD_model_regr$save_object(file=here("fits","BD_model_regr.RDS"))
#load regression model:
BD_model_regr <- readRDS(here("fits","BD_model_regr.RDS"))
#extract posterior:
post_BD_model <- as_draws_matrix(BD_model_regr)
# trace plots:
trace_plot1 <- mcmc_trace(post_BD_model, pars=c("y0_initial[1]","k1","k2"))
trace_plot2 <- mcmc_trace(post_BD_model, pars=c("sigma[1]","sigma[2]","sigma[3]"))
trace_plot3 <- mcmc_trace(post_BD_model, pars=c("Omega[1,2]","Omega[1,3]","Omega[2,3]"))
trace_plot1/trace_plot2/trace_plot3
# diagnosis via shinystan:
shinystan::launch_shinystan(BD_model_regr)
# Comparison of priors and posteriors
samp_BD <- BD_model_regr$draws(format="draws_df")
post_chol_BD <- samp_BD %>% 
  select(!ends_with('prior')) %>% 
  select(!starts_with(".")) %>% 
  select(-"lp__") %>% 
  select(!contains("L_p[")) 

prior_chol_BD <- samp_BD %>%
  select(ends_with("prior")) %>%
  rename_with(~ gsub("_prior", "", .x, fixed = TRUE))

k1_plot <- ggplot(data=prior_chol_BD, aes(x=k1))+
  geom_density(fill="lightblue")+
  geom_density(data=post_chol_BD, aes(x=k1), fill="turquoise")+
  labs(x=expression(k[1]), y="density")+
  xlim(c(0,1))
k2_plot <- ggplot(data=prior_chol_BD, aes(x=k2))+
  geom_density(fill="lightblue")+
  geom_density(data=post_chol_BD, aes(x=k2), fill="turquoise")+
  labs(x=expression(k[2]), y="")+  
  xlim(c(0,2))
A0_plot <-  ggplot(data=prior_chol_BD, aes(x=A0))+
  geom_density(fill="lightblue")+
  geom_density(data=samp_BD, aes(x=`y0_initial[1]`), fill="turquoise")+
  labs(x=expression(A[0]), y="")+
  xlim(c(0,2))

k1_plot+k2_plot + A0_plot
#same for sigma:
sigma1_plot <- ggplot(data = prior_chol_BD, aes(x=sigma))+
  geom_density(fill='lightblue')+
  geom_density(data=post_chol_BD, aes(x=`sigma[1]`), fill='turquoise')+
  xlim(c(0,0.5))+
  labs(x=expression(sigma[A]), y="")

sigma2_plot <- ggplot(data = prior_chol_BD, aes(x=sigma))+
  geom_density(fill='lightblue')+
  geom_density(data=post_chol_BD, aes(x=`sigma[2]`), fill='turquoise')+
  xlim(c(0,0.5))+
  labs(x=expression(sigma[B]),y="")

#Same for correlation coefficients experimental variance-covariance matrix:
LKJ_2 <- 
  rethinking::rlkjcorr(1e6, K = 3, eta = 2) %>%
  data.frame()

Omega_AB <- LKJ_2 %>% 
  ggplot(aes(x = X2)) +
  geom_density(color = "black", fill = "lightblue", alpha = 2/3, adjust = 1/2)+
  geom_density(data=post_chol_BD, aes(x=cor_A_B), fill="turquoise", alpha=0.7)+
  labs(x=expression(rho[AB]), y="density")

Omega_AC <- LKJ_2 %>% 
  ggplot(aes(x = X2)) +
  geom_density(color = "black", fill = "lightblue", alpha = 2/3, adjust = 1/2)+
  geom_density(data=post_chol_BD, aes(x=cor_A_C), fill="turquoise", alpha=0.7)+
  labs(x=expression(rho[AC]), y="")

Omega_BC <- LKJ_2 %>% 
  ggplot(aes(x = X2)) +
  geom_density(color = "black", fill = "lightblue", alpha = 2/3, adjust = 1/2)+
  geom_density(data=post_chol_BD, aes(x=cor_B_C), fill="turquoise", alpha=0.7)+
  labs(x=expression(rho[BC]), y="")

Omega_AB + Omega_AC + Omega_BC
sigma3_plot <- ggplot(data = prior_chol_BD, aes(x=sigma))+
  geom_density(fill='lightblue')+
  geom_density(data=post_chol_BD, aes(x=`sigma[3]`), fill='turquoise')+
  xlim(c(0,0.5))+
  labs(x=expression(sigma[C]),y="")

sigma1_plot+sigma2_plot+sigma3_plot

# use bayesplot to display posteriors:
areas1 <- mcmc_areas(post_BD_model, pars = c("k1","k2", "y0[1]"))+ggplot2::scale_y_discrete(labels=c(expression(paste("k"[1])),expression(paste("k"[2])),expression(paste("A"[0]))))+labs(x="parameter estimates (arbitrary values)")
areas2 <- mcmc_areas(post_BD_model, pars = c("sigma[1]","sigma[2]","sigma[3]")) + ggplot2::scale_y_discrete(labels=c(expression(sigma[A]),expression(sigma[B]),expression(sigma[C])))+labs(x="experimental standard deviations estimates (arbitrary values)")
areas1/areas2
# use GGally to plot pairs and correlation plot:
post_BD_model <- posterior::as_draws_matrix(BD_model_regr) %>% as_tibble() %>% dplyr::select(-lp__) %>% dplyr::select(!starts_with("L_p")) %>% dplyr::select(!starts_with("y_hat")) %>% dplyr::select(!starts_with("y0[")) %>% dplyr::select(-"Omega[1,1]") %>% dplyr::select(-"Omega[2,2]") %>% dplyr::select(-"Omega[3,3]") %>% dplyr::select(-"Omega[1,2]") %>% dplyr::select(-"Omega[1,3]") %>% dplyr::select(-"Omega[2,3]")
cor_BD <- dplyr::select(post_BD,k1:"Omega[3,2]")
# change the names of the columns to be displayed in the panels
cor_BD <- set_names(cor_BD, c(expression(paste("k"[1])), 
                              expression(paste("k"[2])), expression(paste("A"[0])), expression(sigma[A]), expression(sigma[B]), expression(sigma[C])), expression(rho[AB]),expression(rho[AC]),expression(rho[BC]))

# use ggally for a pairs plot
(cor_BD_plot <- cor_BD  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 2, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 9, color = "red"),
          strip.text.y = element_text(size = 9, color = "red"))+
    theme(axis.text.x = element_text(angle=70, hjust=1))
) 
# prepare residuals and QQ plots:
t_pred = seq(from =min(BD_df$time), to = max(BD_df$time), length.out=40)
y_hat <- BD_model_regr %>% 
  spread_draws(c(y_hat)[time_index,component]) %>% 
  mean_qi() %>% 
  mutate(time=time[time_index]) # connect t_pred to time_index
BD_obs_long <- pivot_longer(BD_df, cols = c("A", "B", "C"), names_to = "species")
y_hat$y_measured <- BD_obs_long$value # add the measured values to the predicted values so that resids can be calculated
y_hat$species <- BD_obs_long$species
y_hat <- y_hat %>% mutate(resid=y_measured-y_hat)
# make the plot:
resid_plot_BD <- y_hat %>% ggplot(aes(x=time, y=resid))+
  geom_point()+
  facet_wrap(~species)+
  geom_segment(aes(xend=time, yend=0), color="blue", alpha=0.2)+
  geom_hline(yintercept = 0, lty=2)+
  labs(x="time (min)", y="residuals")
# make the QQ plot:
qqplot_BD <- y_hat %>% ggplot(aes(sample=resid))+
  stat_qq() +
  stat_qq_line()+
  facet_wrap(~species)

resid_plot_BD/qqplot_BD

# calculation of values to compare AthenaVisual Plus with Stan
# results copied from AthenaVisual Plus output
SigmaAA <- 0.01265393/16
SigmaBB <- 0.02979971/16
SigmaCC <- 0.01253283/16
SigmaAB <- -0.008172781/16
SigmaAC <- 0.00594/16
SigmaBC <- 0.005541477/16
sigma_A <-  round(sqrt(SigmaAA),3)
sigma_B <-  round(sqrt(SigmaBB),3)
sigma_C <- round(sqrt(SigmaCC),3)
# results from Stan:
BD_summary <- BD_model_regr$summary()
sigma_A_Stan <- BD_summary[5,2]
sigma_B_Stan <- BD_summary[6,2]
sigma_C_Stan <- BD_summary[7,2]

# correlation matrix for experimental variance-covariance matrix from Athena:
rhoAB <- SigmaAB/(sqrt(SigmaAA)*sqrt(SigmaBB))
rhoAC <- SigmaAC/(sqrt(SigmaAA)*sqrt(SigmaCC))
rhoBC <- SigmaBC/(sqrt(SigmaBB)*sqrt(SigmaCC))
rho_matrix <- matrix(c(1,rhoAB,rhoAC,rhoAB,1,rhoBC,rhoAC,rhoBC,1),3,3, byrow = TRUE)
rho_matrix

#Plot of the correlation coefficients of the experimental variance-covariance matrix using bayesplot: 
mcmc_plot1 <- mcmc_hist(post_BD_model, pars = c("Omega[2,1]"))+geom_vline(xintercept = -0.42, color="red", linewidth=1.5)+labs(x=expression(rho[AB]), y="density")
mcmc_plot2 <- mcmc_hist(post_BD_model, pars = c("Omega[3,1]"))+geom_vline(xintercept = 0.47, color="red", linewidth=1.5)+ labs(x=expression(rho[AC]))
mcmc_plot3 <- mcmc_hist(post_BD_model, pars = c("Omega[3,2]"))+geom_vline(xintercept = 0.29, color="red", linewidth=1.5)+ labs(x=expression(rho[BC]))
mcmc_plot1 + mcmc_plot2 + mcmc_plot3

# nonlinear least-squares regression:
# ODE model function
ABC_func <- function(parms, time) {
  # parameters passed through a list called parms
  A0 = parms[["A0"]]
  k1 = parms[["k1"]]
  k2 = parms[["k2"]]
  
  cA = A0 * exp(-k1 * time)
  cB = (k1 * A0 / (k2 - k1)) * (exp(-k1 * time)-exp(-k2 * time))
  cC = A0 * (1 + ((k1 * exp(-k2 * time) - k2 * exp(-k1 * time)) / (k2 - k1)))
# return function evaluations as vector
  return(c(cA, cB, cC))
}
#starting parameter values for the regression routine:
parms = c(A0 = 1, k1 = 0.08, k2 = 0.3)
# observation vector
expconc = BD_obs$value
# regression routine from gslnls:
fitval_BD <- gsl_nls(
  fn = ABC_func,
  y = expconc,
  start = parms,
  algorithm = "lmaccel",
  time = BD_df$time 
)
# calculation of the parameter confidence intervals
coef_nls_CI <- confint(fitval_BD, level = 0.95, method = "asymptotic")
# comparison of parameter intervals:
#results from Stan:
coef_stan_df <- BD_model_regr$draws(format="df") 
q95_k1 <- quantile(coef_stan_df$k1,probs=c(0.025,0.975))
q95_k2 <- quantile(coef_stan_df$k2,probs=c(0.025,0.975))
q95_c0 <- quantile(coef_stan_df$`y0_initial[1]`,probs=c(0.025,0.975))
coef_stan <- tibble(term=c("A0","k1","k2"), 
            estimate=c(mean(coef_stan_df$`y0_initial[1]`),mean(coef_stan_df$k1),mean(coef_stan_df$k2)),
            CI_low=c(q95_c0[1],q95_k1[1],q95_k2[1]),
            CI_high=c(q95_c0[2],q95_k1[2],q95_k2[2]),
            method=c("Stan", "Stan", "Stan"))
# results from Athena
coef_athena <- tibble(term=c("A0","k1","k2"),
                      estimate=c(0.9949783,0.2075739,0.4962369),
                      CI_low=c((0.9949783-0.01524),(0.2075739-0.01174),(0.4962-0.06928)),
                      CI_high=c((0.9949783+0.01524),(0.2075739+0.01174),(0.4962+0.06928)),
                      method=c("Athena","Athena","Athena"))
# results from LS:
coef_nls <- tibble(term=c("A0","k1","k2"),
                   estimate=c(coef(summary(fitval_BD))[1,1],coef(summary(fitval_BD))[2,1],coef(summary(fitval_BD))[3,1]),
                   CI_low=c(coef_nls_CI[1,1],coef_nls_CI[2,1],coef_nls_CI[3,1]),
                   CI_high=c(coef_nls_CI[1,2],coef_nls_CI[2,2],coef_nls_CI[3,2]),
                   method=c("nls","nls","nls"))
coef_bind <- bind_rows(coef_stan, coef_athena, coef_nls)
ggplot(coef_bind, aes(method, estimate, col=method))+
  geom_point(size=3)+
  geom_linerange(aes(ymin=CI_low, ymax=CI_high))+
  facet_wrap(~term, scales='free_y')+
  theme(legend.position = 'none')
# calculation of fits, confidence and prediction intervals, first Stan result:
t_pred = seq(from =0, to = max(BD_df$time), length.out=40)
y_draws2 <- BD_model_regr %>% 
  spread_draws(c(y_hat2,y_pred)[time_index,component]) %>% 
  mean_qi() %>% 
  mutate(time=t_pred[time_index]) # connect t_pred to time_index
y_draws2_plot <- y_draws2 %>% filter(component==1) %>% ggplot(aes(x=time, y=y_pred))+
  geom_line(color="red", linewidth=1)+
  geom_ribbon(aes(ymin=y_hat2.lower, ymax=y_hat2.upper), fill="red", alpha=0.2) +
  geom_ribbon(aes(ymin=y_pred.lower, ymax=y_pred.upper), fill="red", alpha=0.2)+
  geom_point(data=BD_df, aes(x=time, y=A),color="red", shape=20, size=1)+
  geom_line(data=y_draws2 %>% filter(component==2), aes(x=time, y=y_pred), color="darkgreen", linewidth=1)+
  geom_ribbon(data=y_draws2 %>% filter(component==2),aes(x=time,ymin=y_pred.lower, ymax=y_pred.upper), fill="darkgreen", alpha=0.3)+
  geom_point(data=BD_df, aes(x=time, y=B), color="darkgreen", shape = 17, size=1)+
  geom_line(data=y_draws2 %>% filter(component==3),aes(x=time, y=y_pred), color="blue", linewidth=1)+
  geom_point(data=BD_df, aes(x=time, y=C), shape=15, size=1, color="blue")+
  geom_ribbon(data=y_draws2 %>% filter(component==3),aes(x=time,ymin=y_pred.lower, ymax=y_pred.upper), fill="lightblue", alpha=0.5)+
  labs(x = "Time (arbitrary values)", y = "concentration (arbitrary values)", subtitle = "Bayesian")
#similarly for least-squares:
tpred_ls <- seq(from = 0, to = max(BD_df$time), length.out=50)
BD_pred_ls <- as.data.frame(predict(fitval_BD, newdata = list(time = tpred_ls), interval = "prediction", level = 0.95))
BD_pred_ls <- cbind(time = rep(tpred_ls, times = 3), species = rep(c("A", "B", "C"), each = 50), BD_pred_ls)
y_ls_plot <- ggplot(BD_pred_ls, aes(x = time, color = species)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species), color = NA, alpha = 0.3) +
  geom_line(aes(y = fit)) +
  geom_point(data = BD_obs, aes(y = value, shape = species)) +
  labs(x = "Time (arbitrary values)", y = "", subtitle = "nls")

y_draws2_plot + y_ls_plot

# Pinene data and analysis

pinene_data1 <- read.csv(file=here("data","pinene1.csv"), header=TRUE, sep=",")
# the next line puts the order such that the same time is listed for all components
pinene_data1_long<- pinene_data1 %>% pivot_longer(cols=c("pinene","dipentene","alloocimene","pyronene","dimer"), names_to = "species") 
# the next line with cols_vary = slowest puts the order such that times are listed per component
pinene_data1_long2 <- pinene_data1 %>% pivot_longer(cols=c("pinene","dipentene","alloocimene","pyronene","dimer"), names_to = "species", cols_vary = "slowest")

#Plot of the pinene data:
(pinene_plot <- pinene_data1_long2 %>% ggplot(aes(x=time, y=value, shape=species,colour=species))+
    geom_point(size=2)+
    labs(x="time (days)", y="mol %")
)
# compile the model:
pinene_model <- cmdstan_model(here("fits","pinene.stan"))
# prepare for sampling:
t_pred = seq(from =0.001, to = max(pinene_data1$time), length.out=40)
inits_chain1 <- list(k1=0.08,k2=0.04,k3=0.05,k4=0.03, k5=0.03)
inits_chain2 <- list(k1=0.08,k2=0.04,k3=0.05,k4=0.03, k5=0.03)
inits_chain3 <- list(k1=0.08,k2=0.04,k3=0.05,k4=0.03, k5=0.03)
inits_chain4 <- list(k1=0.08,k2=0.04,k3=0.05,k4=0.03, k5=0.03)

pinene_data_stan <- list(
  N_t = length(pinene_data1$time)-1,
  N_r=5, # number of reactions
  N_m=5, # number of measured components
  ts = pinene_data1$time[-1],
  y = pinene_data1[-1,2:6],
  t0=0,
  y0=100,
  t_pred = t_pred,
  N_p=length(t_pred),
  lkj_df = 2
)

# sampling from the compiled model:
pinene_regr <- pinene_model$sample(
  data = pinene_data_stan,
  seed = 42L,
  refresh = 0,                # silence progress notification
  parallel_chains=4,
  iter_warmup = 4000,
  iter_sampling = 8000,
  #  max_treedepth = 15,
  #  adapt_delta = 0.95,
  init = list(inits_chain1,inits_chain2,inits_chain3,inits_chain4)
)
# store the regression result:
pinene_regr$save_object(file=here("fits","pinene_regr.RDS"))
#load the regression result:
pinene_regr <- readRDS(here("fits","pinene_regr.RDS"))
#Bayesplot for displaying posteriors:
pinene_post <- as_draws_df(pinene_regr)
areas_plot1 <- mcmc_areas(pinene_regr$draws(c("k1","k2","k3")))+ggplot2::scale_y_discrete(labels=c(expression(paste("k"[1])),expression(paste("k"[2])),expression(paste("k"[3]))))+labs(x="rate constants value")
areas_plot1a <- mcmc_areas(pinene_regr$draws(c("k4","k5")))+ggplot2::scale_y_discrete(labels=c(expression(paste("k"[4])),expression(paste("k"[5]))))+
  labs(x="rate constants value")
areas_plot2 <- mcmc_areas(pinene_regr$draws(c("sigma[1]","sigma[2]","sigma[3]","sigma[4]","sigma[5]")))+ggplot2::scale_y_discrete(labels=c(expression(sigma[A]),expression(sigma[B]),expression(sigma[C]),expression(sigma[D]),expression(sigma[E])))+labs(x="experimental standard deviations")
areas_plot3 <- mcmc_areas(pinene_regr$draws(c("y0_initial[1]")))+ggplot2::scale_y_discrete(labels=c(expression(paste("A"[0]))))+labs(x="initial value pinene")
(areas_plot1+areas_plot1a)/(areas_plot3+areas_plot2)

t_pred = seq(from =0.001, to = max(pinene_data1$time), length.out=40)

y_draws_pinene<- pinene_regr %>% 
  spread_draws(c(y_pred)[time_index,component]) %>% 
  mean_qi() %>% 
  mutate(time=t_pred[time_index]) # connect t_pred to time_index
#Plot the fit, residuals and QQ plots:
y_draws_pinene_plot <- y_draws_pinene %>% filter(component==1) %>% ggplot(aes(x=time, y=y_pred))+
  geom_line(color="blue", linewidth=1, alpha=0.3)+
  geom_point(data=pinene_data1_long, aes(x=time, y=value,color=species, shape=species))+
  geom_line(data=y_draws_pinene %>% filter(component==2),aes(x=time, y=y_pred), color="darkgreen", linewidth=1)+  
  geom_line(data=y_draws_pinene %>% filter(component==3),aes(x=time, y=y_pred), color="chocolate1", linewidth=1)+
  geom_line(data=y_draws_pinene %>% filter(component==5),aes(x=time, y=y_pred), color="darkgoldenrod1", linewidth=1)+
  geom_line(data=y_draws_pinene %>% filter(component==4),aes(x=time, y=y_pred), color="darkorchid1", linewidth=1)+
  labs(x="time (days)", y="mol%")
#remove the t0 values
ts = pinene_data1$time[-1]

y_hat<- pinene_regr %>% 
  spread_draws(c(y_hat)[time_index,component]) %>% 
  mean_qi() %>% 
  mutate(time=ts[time_index])

y_hat$y_measured <- pinene_data1_long_no_t0$value 
y_hat <- y_hat %>% mutate(resid=y_measured-y_hat)
species <- c("pinene","dipentene","alloocimene","pyronene","dimer")
species_df <- rep(species,8)
y_hat <- cbind(y_hat,species_df)
resids_plot_pinene <- y_hat %>% ggplot(aes(x=time, y=resid))+
  geom_point()+
  facet_wrap(~factor(species_df,c("pinene","dipentene","alloocimene","pyronene","dimer")))+
  geom_segment(aes(xend=time, yend=0), color="blue", alpha=0.2)+
  geom_hline(yintercept = 0, lty=2)+
  labs(x="time (days)", y="residuals")
qqplot_pinene <- y_hat %>% ggplot(aes(sample=resid))+
  stat_qq() +
  stat_qq_line()+
  facet_wrap(~factor(species_df,c("pinene","dipentene","alloocimene","pyronene","dimer")))
y_draws_pinene_plot/(resids_plot_pinene+qqplot_pinene)

pinene_post_corr <- pinene_post %>% select(starts_with("cor"))
# Plot of the correlation coefficients of the experimental variance-covariance matrix:
# use as_labeller() to display greek characters:
my_labeller_two <- as_labeller(
  x = c('cor_A_B' = 'rho[AB]','cor_A_C' = 'rho[AC]','cor_A_D' = 'rho[AD]','cor_A_E' = 'rho[AE]','cor_B_C' = 'rho[BC]','cor_B_D' = 'rho[BD]','cor_B_E' = 'rho[BE]','cor_C_D' = 'rho[CD]','cor_C_E' = 'rho[CE]','cor_D_E' = 'rho[DE]'), 
  default = label_parsed
)
mcmc_omega <- mcmc_hist(pinene_post_corr, facet_args = list(labeller = my_labeller_two))
mcmc_omega

# Parameter correlation for pinene:
pinene_regr <- readRDS(here("fits","pinene_regr.RDS"))

pinene_post <- as_draws_df(pinene_regr)

cor_pinene <- pinene_post %>% select(k1:`sigma[5]`) 
cor_pinene2 <- pinene_post %>% select(starts_with('cor'))
cor_pinene <- cbind(cor_pinene,cor_pinene2) 
cor_coef <- cor(cor_pinene)
colnames(cor_coef) <- c(":k[1]",":k[2]",":k[3]",":k[4]",":k[5]",":A[0]",":sigma[A]",":sigma[B]",":sigma[C]",":sigma[D]",":sigma[E]",":rho[AB]",":rho[AC]",":rho[AD]",":rho[AE]",":rho[BC]",":rho[BD]",":rho[BE]",":rho[CD]",":rho[CE]",":rho[DE]")
rownames(cor_coef) <- c(":k[1]",":k[2]",":k[3]",":k[4]",":k[5]",":A[0]",":sigma[A]",":sigma[B]",":sigma[C]",":sigma[D]",":sigma[E]",":rho[AB]",":rho[AC]",":rho[AD]",":rho[AE]",":rho[BC]",":rho[BD]",":rho[BE]",":rho[CD]",":rho[CE]",":rho[DE]")
corrplot::corrplot(cor(cor_coef))

# Prior-posterior comparison for pinene:
samp_pinene <- pinene_regr$draws(format="draws_df")
post_chol_pinene <- samp_pinene %>% 
  select(!ends_with('prior')) %>% 
  select(!starts_with(".")) %>% 
  select(-"lp__") %>% 
  select(!contains("L_p[")) 

prior_chol_pinene <- samp_pinene %>%
  select(ends_with("prior")) %>%
  rename_with(~ gsub("_prior", "", .x, fixed = TRUE))

k1_plot <- ggplot(data=prior_chol_pinene, aes(x=k1))+
  geom_density(fill="red")+
  geom_density(data=post_chol_pinene, aes(x=k1), fill="turquoise")+
  labs(x=expression(k[1]), y="density")+
  xlim(c(0,0.1))+
  theme(axis.text.x = element_text(angle=70, hjust=1))
k2_plot <- ggplot(data=prior_chol_pinene, aes(x=k2))+
  geom_density(fill="red")+
  geom_density(data=post_chol_pinene, aes(x=k2), fill="turquoise")+
  labs(x=expression(k[2]), y="")+  
  xlim(c(0,0.1))+
  theme(axis.text.x = element_text(angle=70, hjust=1))
k3_plot <- ggplot(data=prior_chol_pinene, aes(x=k3))+
  geom_density(fill="red")+
  geom_density(data=post_chol_pinene, aes(x=k2), fill="turquoise")+
  labs(x=expression(k[3]), y="")+  
  xlim(c(0,0.1))+
  theme(axis.text.x = element_text(angle=70, hjust=1))
k4_plot <- ggplot(data=prior_chol_pinene, aes(x=k4))+
  geom_density(fill="red")+
  geom_density(data=post_chol_pinene, aes(x=k2), fill="turquoise")+
  labs(x=expression(k[4]), y="")+  
  xlim(c(0,0.1))+
  theme(axis.text.x = element_text(angle=70, hjust=1))
k5_plot <- ggplot(data=prior_chol_pinene, aes(x=k5))+
  geom_density(fill="red")+
  geom_density(data=post_chol_pinene, aes(x=k2), fill="turquoise")+
  labs(x=expression(k[5]), y="")+  
  xlim(c(0,0.1))+
  theme(axis.text.x = element_text(angle=70, hjust=1))
y0_plot <-  ggplot(data=prior_chol_pinene, aes(x=`y0`))+
  geom_density(fill="red")+
  geom_density(data=samp_pinene, aes(x=`y0_initial[1]`), fill="turquoise")+
  labs(x=expression(y[0]), y="")+
  xlim(c(0,200))+
  theme(axis.text.x = element_text(angle=70, hjust=1))

(k1_plot+k2_plot + k3_plot)/(k4_plot+k5_plot+y0_plot)

#Least-squares regression of the pinene data set:
pinene_ls <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dy1 <-  -exp(log(k1)) * y1 - exp(log(k2))*y1
    dy2 <-  exp(log(k1)) * y1
    dy3 <- exp(log(k2))*y1 - exp(log(k3))*y3 -2*exp(log(k4))*y3*y3 + 2*exp(log(k5))*y5
    dy4 <- exp(log(k3))*y3
    dy5 <- exp(log(k4))*y3*y3 - exp(log(k5))*y5
    list(c(dy1, dy2,dy3, dy4,dy5))
  })
}

# function that calls ode integrator and returns concentration values
ode_func <- function(parms, time) {
predconc <- deSolve::ode(y = c(parms["y1"],y2=0.0,y3=0.0,y4=0.0,y5=0.0), times = time, parms = parms[c("k1","k2","k3","k4","k5")], func = pinene_ls)
# return concentrations as vector, this returns the values for y1,y2, y3 for each time point
  return(c(predconc[, c("y1", "y2", "y3","y4","y5")])) }
# observation vector needs this format:
expconc_pinene = pinene_data1_long2$value
t_obs <- pinene_data1$time
# starting values for regression (with y0[1] as parameter, the other starting values are put at zero):
y0 <- c(y1 = 100.0)
parms = c(y0, k1 = 0.08, k2 = 0.04, k3=0.03, k4=0.04, k5=0.06)
fit_pinene_ls <- gsl_nls(
  fn = ode_func,
  y = expconc_pinene,
  start = parms,
  algorithm = "lmaccel",
  time = t_obs
)
#Summary of least-squares of pinene regression
summary(fit_pinene_ls)
# 95% parameter confidence levels from least-squares:
confint(fit_pinene_ls,method = "asymptotic", level = 0.95)
# plot fitted trends + prediction intervals resulting from least-squares regression
tpred <- seq(from = 0, to = 25, length = 50)
pinene_pred <- as.data.frame(predict(fit_pinene_ls, newdata = list(time = tpred), interval = "prediction", level = 0.95))
pinene_pred <- cbind(time = rep(tpred, times = 5), species = rep(c("pinene", "dipentene", "alloocymene","pinonene","dimer"), each = 50), pinene_pred)

ggplot(pinene_pred, aes(x = time, y=fit)) +
  geom_point(data = pinene_data1_long, aes(y = value, color=species)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species),alpha = 0.3) +
  labs(x = "time (days)", y = "mol %")+
  guides(fill = FALSE, size = FALSE)

# Comparison of confidence intervals from Stan and nls:
coef_nls_CI_pinene <- confint(fit_pinene_ls, level = 0.95, method = "asymptotic")
# results from Stan:
coef_stan_pinene <- pinene_regr$draws(format="df") 
q95_k1 <- quantile(coef_stan_pinene$k1,probs=c(0.025,0.975))
q95_k2 <- quantile(coef_stan_pinene$k2,probs=c(0.025,0.975))
q95_k3 <- quantile(coef_stan_pinene$k3,probs=c(0.025,0.975))
q95_k4 <- quantile(coef_stan_pinene$k4,probs=c(0.025,0.975))
q95_k5 <- quantile(coef_stan_pinene$k5,probs=c(0.025,0.975))
q95_c0 <- quantile(coef_stan_pinene$`y0_initial[1]`,probs=c(0.025,0.975))
coef_stan <- tibble(term=c("A0","k1","k2","k3","k4","k5"), 
                    estimate=c(mean(coef_stan_pinene$`y0_initial[1]`),mean(coef_stan_pinene$k1),mean(coef_stan_pinene$k2),mean(coef_stan_pinene$k3),mean(coef_stan_pinene$k4),mean(coef_stan_pinene$k5)),
                    CI_low=c(q95_c0[1],q95_k1[1],q95_k2[1],q95_k3[1],q95_k4[1],q95_k5[1]),
                    CI_high=c(q95_c0[2],q95_k1[2],q95_k2[2],q95_k3[2],q95_k4[2],q95_k5[2]),
                    method=c("Stan", "Stan", "Stan","Stan","Stan","Stan"))
# results from LS:
coef_nls <- tibble(term=c("A0","k1","k2","k3","k4","k5"),
                   estimate=c(coef(summary(fit_pinene_ls))[1,1],coef(summary(fit_pinene_ls))[2,1],coef(summary(fit_pinene_ls))[3,1],coef(summary(fit_pinene_ls))[4,1],coef(summary(fit_pinene_ls))[5,1],coef(summary(fit_pinene_ls))[6,1]),
                   CI_low=c(coef_nls_CI_pinene[1,1],coef_nls_CI_pinene[2,1],coef_nls_CI_pinene[3,1],coef_nls_CI_pinene[4,1],coef_nls_CI_pinene[5,1],coef_nls_CI_pinene[6,1]),
                   CI_high=c(coef_nls_CI_pinene[1,2],coef_nls_CI_pinene[2,2],coef_nls_CI_pinene[3,2],coef_nls_CI_pinene[4,2],coef_nls_CI_pinene[5,2],coef_nls_CI_pinene[6,2]),
                   method=c("nls","nls","nls","nls","nls","nls"))

coef_bind <- bind_rows(coef_stan, coef_nls)

ggplot(coef_bind, aes(method, estimate, col=method))+
  geom_point(size=3)+
  geom_linerange(aes(ymin=CI_low, ymax=CI_high))+
  facet_wrap(~term, scales='free_y')+
  theme(legend.position = 'none')

t_pred = seq(from =0.001, to = max(pinene_data1$time), length.out=40)
# y_hat2 contains the variation around the mean (variation of the parameters), y_pred the variation including the experimental standard deviations, both 95% levels:
y_draws_pinene<- pinene_regr %>% 
  spread_draws(c(y_hat2,y_pred)[time_index,component]) %>% 
  mean_qi() %>% 
  mutate(time=t_pred[time_index]) # connect t_pred to time_index

#calculations for the least-squares regression:
tpred_ls <- seq(from = 0, to = 25, length = 50)

# ls calculation for 95% prediction:
pinene_pred <- as.data.frame(predict(fit_pinene_ls, newdata = list(time = tpred_ls), interval = "prediction", level = 0.95))
pinene_pred <- cbind(time = rep(tpred_ls, times = 5), species = rep(c("pinene", "dipentene", "alloocimene","pyronene","dimer"), each = 50), pinene_pred)
# ls calculation for 95% confidence:
pinene_conf <- as.data.frame(predict(fit_pinene_ls, newdata = list(time = tpred_ls), interval = "confidence", level = 0.95))
pinene_conf <- cbind(time = rep(tpred_ls, times = 5), species = rep(c("pinene", "dipentene", "alloocimene","pyronene","dimer"), each = 50), pinene_conf)
# Comparison of Stan and least-squares result:
confint_pinene <- ggplot(data=pinene_conf %>% filter(species=="pinene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha = 0.5)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==1), aes(y=y_hat2,ymin=y_hat2.lower, ymax=y_hat2.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=pinene))+
  labs(x="time (days)", y="mol% pinene")
confint_dipentene <- ggplot(data=pinene_conf %>% filter(species=="dipentene"), aes(x=time,y=fit))+
  geom_line()+  
  geom_ribbon(aes(y=fit, ymin = lwr, ymax = upr),fill="blue", alpha = 0.5)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==2), aes(y=y_hat2,ymin=y_hat2.lower, ymax=y_hat2.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=dipentene))+
  labs(x="time (days)", y="mol% dipentene")
confint_alloo <- ggplot(data=pinene_conf %>% filter(species=="alloocimene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha=0.5)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==3), aes(y=y_hat2,ymin=y_hat2.lower, ymax=y_hat2.upper), fill="orange", alpha=0.3)+
  geom_point(data=pinene_data1,aes(y=alloocimene))+
  labs(x="time (days)", y="mol% alloocimene")
confint_pyro <- ggplot(data=pinene_conf %>% filter(species=="pyronene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(y=fit, ymin = lwr, ymax = upr),fill="blue", alpha = 0.4)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==4), aes(y=y_hat2,ymin=y_hat2.lower, ymax=y_hat2.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=pyronene))+
  labs(x="time (days)", y="mol% pyronene")
confint_dimer <- ggplot(data=pinene_conf %>% filter(species=="dimer"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(y=fit, ymin = lwr, ymax = upr),fill="blue", alpha = 0.4)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==5), aes(y=y_hat2,ymin=y_hat2.lower, ymax=y_hat2.upper), fill="orange", alpha=0.7)+
  geom_point(data=pinene_data1,aes(y=dimer))+
  labs(x="time (days)", y="mol% dimer")

(confint_pinene+confint_dipentene+confint_alloo+confint_pyro+confint_dimer)

#Same for prediction intervals:
pred_pinene <- ggplot(data=pinene_pred %>% filter(species=="pinene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha = 0.7)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==1), aes(y=y_hat2,ymin=y_pred.lower, ymax=y_pred.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=pinene))+
  labs(x="time (days)", y="mol% pinene")

pred_dipentene <- ggplot(data=pinene_pred %>% filter(species=="dipentene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha = 0.7)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==2), aes(y=y_hat2,ymin=y_pred.lower, ymax=y_pred.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=dipentene))+
  labs(x="time (days)", y="mol% dipentene")

pred_alloo <- ggplot(data=pinene_pred %>% filter(species=="alloocimene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha = 0.7)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==3), aes(y=y_hat2,ymin=y_pred.lower, ymax=y_pred.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=alloocimene))+
  labs(x="time (days)", y="mol% alloocimene")

pred_pyro <- ggplot(data=pinene_pred %>% filter(species=="pyronene"), aes(x=time, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha = 0.7)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==4), aes(y=y_hat2,ymin=y_pred.lower, ymax=y_pred.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=pyronene))+
  labs(x="time (days)", y="mol% pyronene")

pred_dimer <- ggplot(data=pinene_pred %>% filter(species=="dimer"), aes(x=time, y=fit))+
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill="blue", alpha = 0.7)+
  geom_ribbon(data=y_draws_pinene %>% filter(component==5), aes(y=y_hat2,ymin=y_pred.lower, ymax=y_pred.upper), fill="orange", alpha=0.4)+
  geom_point(data=pinene_data1,aes(y=dimer))+
  labs(x="time (days)", y="mol% dimer")

(pred_pinene+pred_dipentene+pred_alloo+pred_pyro+pred_dimer)
