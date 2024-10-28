expand_conditions<-function(Js_vec,Ns_vec,psi_vec,alpha_vec){
  #J: Number of items
  #N: Sample size
  #psi: SD of LV distribution 
  #alpha: Skewness parameter
  expand.grid(J = Js_vec, N = Ns_vec, psi = psi_vec, alpha = alpha_vec)
}

get_forms_list<-function(vals_1PL, Js_vec){
  # Function to get items when "sampling down" and creating short forms
  #vals_1PL: parameter estimates from mirt object
  #Js_vec: Vector of form length
  d1_1PL = vals_1PL %>% dplyr::filter(name == "d1") %>% arrange(-value)
  form_list = lapply(Js_vec[-length(Js_vec)], function(J){
    foo = qnorm(seq(0,1, len = J))
    foo[foo == -Inf] = -10; foo[foo==Inf] = 10
    return(sapply(1:J, function(q){d1_1PL$item[which.min(abs(d1_1PL$value-foo[q]))]}))
  })
  form_list[[length(form_list) + 1]] = d1_1PL$item
  names(form_list) = paste0("J",Js_vec)
  return(form_list)
}


simulation_results<-function(conditions_df, my_seed, fit_1PL, form_list){
  
  results<-lapply(1:nrow(conditions_df), function(i){
    
    
    vals_1PL = mirt::mod2values(fit_1PL)
    
    #Extract the simulation parameter conditions
    J_i = conditions_df$J[i]; N_i = conditions_df$N[i]; Psi_i = conditions_df$psi[i]; alpha_i = conditions_df$alpha[i]
    
    # If empirical prior, initiate a new fitting object
    fit_1PL_i = fit_1PL
    if(!is.na(Psi_i) & Psi_i != 1){
      vals_1PL_i = vals_1PL
      vals_1PL_i$est = F
      vals_1PL_i$value[vals_1PL_i$name == "COV_11"] = Psi_i
      fit_1PL_i =  mirt::mirt(data= psy_list %>% purrr::pluck("data") %>% dplyr::select(starts_with("PS")) %>% dplyr::select(-"PS033"), technical = list(NCYCLES = 2000),  TOL = 1E-5, constrain = list(avec_1PL), pars = vals_1PL_i)
    }
    
    #Conduct the simulation
    with_progress({
      p<-progressor(along = 1:Reps)
      ret_i <- future.apply::future_lapply(1:Reps,function(rep){
        
        # Simulate latent scores
        etas_i = sn::rsn(n=N_i, alpha = alpha_i) %>% as.matrix()
        
        #Simulate item responses
        simdat = mirt::simdata(model = fit_1PL_i, Theta = etas_i) %>% data.frame()
        #Shorten the form, as necessary
        vals_1PL_i = mirt::mod2values(fit_1PL_i)
        vals_1PL_i$est = F
        if(J_i!=46){
          items_kept_i = form_list %>% purrr::pluck(paste0("J",J_i))
          items_removed_i = setdiff(form_list[["J46"]], items_kept_i)
          simdat = simdat %>% dplyr::select(-dplyr::any_of(items_removed_i))
          vals_1PL_i = vals_1PL_i%>% dplyr::filter(!(item %in% items_removed_i)) %>%
            dplyr::mutate(parnum = 1:n())
        }
        
        # Obtain empirical hyperparameters, if called for in simulation conditoin
        if(is.na(Psi_i)){
          vals_1PL_i$est[vals_1PL_i$name %in% c("COV_11","MEAN_1")] = T
        }
        
        #Obtain MLE, EAP, and PV scores
        fit_1PL_i =  mirt::mirt(data= simdat, technical = list(NCYCLES = 2000),  TOL = 1E-5, pars = vals_1PL_i)
        MLE_scores = fscores(fit_1PL_i, response.pattern = simdat,method ="ML") %>% data.frame() %>% purrr::pluck("F1")
        EAP_scores = fscores(fit_1PL_i, response.pattern = simdat, method = "EAP") %>% data.frame() %>% purrr::pluck("F1")
        PVs_list = mirt::fscores(fit_1PL_i, response.pattern = simdat, plausible.draws = M)
        
        # Obtain parameter estimates
        mu_MLE = mean(MLE_scores); sd_MLE = sd(MLE_scores); PrOff_MLE = mean(MLE_scores>qnorm(.95)); skew_MLE = moments::skewness(MLE_scores)
        mu_EAP = mean(EAP_scores); sd_EAP = sd(EAP_scores); PrOff_EAP = mean(EAP_scores>qnorm(.95)); skew_EAP = moments::skewness(EAP_scores)
        
        # Rubin's rules (1987): Take mean of posterior predictive
        mu_PV = 0; sd_PV = 0; PrOff_PV = 0; skew_PV = 0
        for(m in 1:M){
          mu_PV = mean(PVs_list[[m]][,1]) + mu_PV
          sd_PV = sd(PVs_list[[m]][,1]) + sd_PV
          PrOff_PV =  mean(PVs_list[[m]][,1]>qnorm(.95)) + PrOff_PV
          skew_PV =moments::skewness(PVs_list[[m]][,1]) + skew_PV
        }
        mu_PV = mu_PV/M; sd_PV = sd_PV/M; PrOff_PV =PrOff_PV/M; skew_PV = skew_PV/M
        
        # Return results
        ret_rep <- data.frame(
          method = c("MLE","EAP","PV"),
          mean_est = c(mu_MLE, mu_EAP, mu_PV),
          sd_est = c(sd_MLE, sd_EAP, sd_PV),
          skew_est = c(skew_MLE, skew_EAP, skew_PV),
          PrOff = c(PrOff_MLE, PrOff_EAP, PrOff_PV)
        ) %>% dplyr::mutate(replication = rep) %>% dplyr::relocate(replication)
        
        p(sprintf("Rep=%g", rep))
        return(ret_rep)
      }, future.seed = my_seed) %>% dplyr::bind_rows()
      
    }) #end with progress
    
    # Summarize simulation results and return summary
    ret_i = ret_i %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(mean_est = mean(mean_est),
                       sd_est = mean(sd_est),
                       skew_est = mean(skew_est),
                       PrOff = mean(PrOff)
      )
    ret_i = ret_i %>% dplyr::mutate(J = J_i, N = N_i, psi = Psi_i, alpha = alpha_i) %>% dplyr::relocate(J:alpha)
    return(ret_i)
    
  }) %>% dplyr::bind_rows()
  
}