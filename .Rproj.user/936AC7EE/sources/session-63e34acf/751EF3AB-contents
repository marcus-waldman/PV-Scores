rm(list = ls())

require(tidyverse)
require(mirt)
require(sn)
require(moments)
require(future)
require(future.apply)
require(progressr)
require(ggthemes)

# Initialize parallel environment
future::plan("multisession", workers = 16)
progressr::handlers("progress")
plan()

# Source in utility functions
source("utils.R")

# Specify file directory for plots
if(!exists("plots")){dir.create("plots")}
if(!exists("results")){dir.create("results")}

# Read in GSED-PF data collected in Nebraska and construct total scores
psy_list = readr::read_rds("data/psy_list.rds")
dat = psy_list$data %>% dplyr::select(starts_with("PS")) %>% dplyr::select(-"PS033")
tot_scores = dat %>% dplyr::select(starts_with("PS")) %>% apply(1,"mean", na.rm = T)
hist(tot_scores) # looks good

# Fit a 1PL (Rasch) Model
vals_1PL =  mirt::mirt(data= psy_list %>% purrr::pluck("data") %>% dplyr::select(starts_with("PS")) %>% dplyr::select(-"PS033"), technical = list(NCYCLES = 2000),  TOL = 1E-5, pars = "values")
vals_1PL$value[vals_1PL$name == "a1"] = 1
avec_1PL = vals_1PL$parnum[vals_1PL$name=="a1"]
fit_1PL =  mirt::mirt(data= psy_list %>% purrr::pluck("data") %>% dplyr::select(starts_with("PS")) %>% dplyr::select(-"PS033"), technical = list(NCYCLES = 2000),  TOL = 1E-5, constrain = list(avec_1PL))
summary(fit_1PL) # looks good


### Plot total scores and give an example of 1:1 correspondence ###
# Impute missing indicator data for demonstration purposes
dat2 = mirt::imputeMissing(fit_1PL, Theta = fscores(fit_1PL)) 
tmp = data.frame(tot_score = apply(dat2, 1, "sum")) %>% 
  dplyr::mutate(eap_score = mirt::fscores(fit_1PL, method = "EAPsum", response.pattern = dat2)[,1])
#1-to-1 correspondence example
plt = ggplot(tmp, aes(x=tot_score, y = eap_score)) + 
  geom_point() + 
  theme_bw() +
  labs(x = "Total score", y = expression("Scaled score, " ~ hat(eta)))
plt
ggsave(plt, filename = "plots/one-to-one-1pl.png", height = 3, width = 3)
# total scores histogram
plt = ggplot(tmp, aes(tot_score)) + geom_histogram() + labs(x = "Total Score", y = NULL) + theme_bw()
ggsave(plt, filename = "plots/total-psy-scores.png", height = 3, width = 3)


### Simulation 1: Ideal conditions: prior dist eq. pop. LV dist: ###

# Conduct simulation
M = 1000 #Draws from posterior [eta | y; delta]
Reps = 2^10 #Simulation replications
# Define simulation conditions
conditions_df = expand_conditions(Js_vec = c(46), Ns_vec = c(500), psi_vec = 1, alpha_vec = 0)
form_list = get_forms_list(vals_1PL = vals_1PL, Js_vec = sort(unique(conditions_df$Js_vec)))
res1<-simulation_results(conditions_df = conditions_df, my_seed = 42, fit_1PL = fit_1PL, form_list = form_list)

# Analyze results
# Wrange to get into format to make available in ggplot2
longtmp = res1 %>% 
  tidyr::pivot_longer(mean_est:PrOff, names_to = "statistic", values_to = "est") %>% 
  dplyr::filter(!is.nan(est), is.finite(est)) %>% #Unbounded likelihood
  dplyr::mutate(param = NA,
                param = ifelse(statistic == "mean_est" | statistic == "skew_est", 0, param), 
                param = ifelse(statistic == "sd_est", 1, param), 
                param = ifelse(statistic == "PrOff", .05, param), 
                stxt = 3,
                stxt = ifelse(statistic=="mean_est", 0, stxt), 
                stxt = ifelse(statistic=="sd_est",1, stxt), 
                stxt = ifelse(statistic=="skew_est", 2, stxt),
                stxt = factor(stxt, level = 0:3, labels = c("Mean", "SD", "Skewness", "Pr(eta>1.65)"), ordered = T),
                estimator = plyr::mapvalues(method, from = c("MLE","EAP","PV"), to = 0:2), 
                estimator = factor(estimator, level = 0:2, labels = c("MLE","EAP","PV")), 
                form = factor(J, levels = c(10,46), labels = c("Short", "Long"), ordered = T)) %>% 
  dplyr::mutate(bias =  est-param)
# Plot in ggplot 2
plt = ggplot(longtmp, aes(x = form, y = bias)) + 
  facet_grid(stxt~estimator) + 
  geom_hline(yintercept = 0) + 
  geom_bar(aes(col = stxt, fill = stxt), stat = "identity", alpha = .5, show.legend = F) + 
  geom_point(shape = 15) +
  theme_minimal() +
  labs(x = "Form Length", y = NULL, title = "Bias (N = 500)", subtitle = "Replications = 1028" ) +
  ggthemes::scale_color_colorblind() + ggthemes::scale_fill_colorblind() + 
  coord_cartesian(ylim = c(-.2,.2))
plt
# Write out the result
ggsave(plt, filename = "plots/Bias-Sim1-N500.png", height = 5, width = 5)



### Simulation 2: Would not informative priors help maleffectss of bayesin shrinkage relieve Bayesian shrinkage: ###

# Conduct simulation
M = 1000 #Draws from posterior [eta | y; delta]
Reps = 2^10 #Simulation replications
# Define simulation conditions
conditions_df = expand_conditions(Js_vec = c(46), Ns_vec = c(500), psi_vec=10, alpha_vec = 0)
form_list = get_forms_list(vals_1PL = vals_1PL, Js_vec = sort(unique(conditions_df$Js_vec)))
res2<-simulation_results(conditions_df = conditions_df, my_seed = 42, fit_1PL = fit_1PL, form_list = form_list)

# Wrangle to get into format to make available in ggplot2
longtmp = res1 %>% 
  dplyr::bind_rows(res2) %>% 
  tidyr::pivot_longer(mean_est:PrOff, names_to = "statistic", values_to = "est") %>% 
  dplyr::filter(!is.nan(est), is.finite(est)) %>% #Unbounded likelihood
  dplyr::mutate(param = NA,
                param = ifelse(statistic == "mean_est" | statistic == "skew_est", 0, param), 
                param = ifelse(statistic == "sd_est", 1, param), 
                param = ifelse(statistic == "PrOff", .05, param), 
                stxt = 3,
                stxt = ifelse(statistic=="mean_est", 0, stxt), 
                stxt = ifelse(statistic=="sd_est",1, stxt), 
                stxt = ifelse(statistic=="skew_est", 2, stxt),
                stxt = factor(stxt, level = 0:3, labels = c("Mean", "SD", "Skewness", "Pr(eta>1.65)"), ordered = T),
                estimator = plyr::mapvalues(method, from = c("MLE","EAP","PV"), to = 0:2), 
                estimator = factor(estimator, level = 0:2, labels = c("MLE","EAP","PV")), 
                form = factor(J, levels = c(10,46), labels = c("Short", "Long"), ordered = T), 
                prior = factor(psi, levels = c(1,10), labels = c("Inform.: N(0,1)", "Noninform.: N(0,10)"), ordered = T)
  )%>% 
  dplyr::mutate(bias =  est-param)
plt = ggplot(longtmp, aes(x = form, y = bias)) + 
  facet_grid(stxt~prior) + 
  geom_hline(yintercept = 0) + 
  geom_bar(aes(col = stxt, fill = stxt), stat = "identity", alpha = .5, show.legend = F) + 
  geom_point(shape = 15) +
  theme_minimal() +
  labs(x = "Form Length", y = NULL, title = "Bias (N = 500)", subtitle = "Replications = 1028" ) +
  scale_color_colorblind() + scale_fill_colorblind() + coord_cartesian(ylim = c(-.2,.2))
plt
# Write out the plots
ggsave(plt, filename = "plots/Bias-Sim2-N500.png", height = 5, width = 5)



# Simulation study 3: What about a skewed distribution? ###

  #Ensure dsn is parameterized as I expect it
    xi = 0
    omega = 1
    alpha = 10; delta = alpha/(sqrt(1+alpha^2))
    
    mu_sn = xi + omega*delta*sqrt(2/pi)
    mean(sn::rsn(1E5,alpha=alpha))#works
    
    var_sn = omega^2 * (1- (2*delta^2/pi) )
    sd_sn = sqrt(var_sn)
    sd(sn::rsn(1E5,alpha=alpha)) #works
    
    num = (delta*sqrt(2/pi))^3
    denom = (1-2*delta^2/pi)^(3/2)
    skew_sn = .5*(4-pi)*(num/denom)
    moments::skewness(sn::rsn(1E5, alpha = alpha)) #works

  # Conduct simulation
    M = 1000 #Draws from posterior [eta | y; delta]
    Reps = 2^10 #Simulation replications
    # Define simulation conditions
    conditions_df = expand_conditions(Js_vec = c(46), Ns_vec = c(500), psi_vec=c(1,10,NA), alpha_vec = 10)
    form_list = get_forms_list(vals_1PL = vals_1PL, Js_vec = sort(unique(conditions_df$Js_vec)))
    res3<-simulation_results(conditions_df = conditions_df, my_seed = 42, fit_1PL = fit_1PL, form_list = form_list)
    
  # Wrangle in preparation for ggplot2
    longtmp = res3 %>% 
      tidyr::pivot_longer(mean_est:PrOff, names_to = "statistic", values_to = "est") %>% 
      dplyr::filter(!is.nan(est), is.finite(est)) %>% 
      dplyr::mutate(param = NA,
                    param = ifelse(statistic == "mean_est", mu_sn, param), 
                    param = ifelse(statistic == "sd_est", sd_sn, param),
                    param =ifelse(statistic == "skew_est", skew_sn, param),
                    param = ifelse(statistic == "PrOff", .1, param), 
                    stxt = 3,
                    stxt = ifelse(statistic=="mean_est", 0, stxt), 
                    stxt = ifelse(statistic=="sd_est",1, stxt), 
                    stxt = ifelse(statistic=="skew_est", 2, stxt),
                    stxt = factor(stxt, level = 0:3, labels = c("Mean", "SD", "Skewness", "Pr(eta>1.65)"), ordered = T),
                    estimator = plyr::mapvalues(method, from = c("MLE","EAP","PV"), to = 0:2), 
                    estimator = factor(estimator, level = 0:2, labels = c("MLE","EAP","PV")), 
                    form = factor(J, levels = c(10,46), labels = c("Short", "Long"), ordered = T), 
                    prior = 0,
                    prior = ifelse(psi==10, 1,prior), 
                    prior = ifelse(is.na(psi), 2, prior), 
                    prior = factor(prior, levels = c(0,1,2), labels = c("Noninform.","Inform.","Empirical"), ordered = T) 
      )%>% 
      dplyr::mutate(bias =  est-param) %>% 
      dplyr::filter(estimator != "MLE")
    plt = ggplot(longtmp, aes(x = form, y = bias)) + 
      facet_grid(stxt~prior+estimator) + 
      geom_hline(yintercept = 0) + 
      geom_bar(aes(col = stxt, fill = stxt), stat = "identity", alpha = .5, show.legend = F) + 
      geom_point(shape = 15) +
      theme_minimal() +
      labs(x = "Form Length", y = NULL, title = "Bias (N = 500)", subtitle = "Replications = 1028" ) +
      scale_color_colorblind() + scale_fill_colorblind() + coord_cartesian(ylim = c(-.2,.2))
    plt
    #Save out the plot
    ggsave(plt, filename = "plots/Bias-Sim3-N500.png", height = 5, width = 8)
    

 
  
