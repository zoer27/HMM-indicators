#HMM Models for 
#Case Study 1
#Ocean temperature from stoplight chart


library(tidyverse)
library(hmmTMB) #for HMM models
library(patchwork) #for plotting



# Data --------------------------------------------------------------------

#reading in stoplight data
stoplight<-read_csv("Data/stoplight2022.csv")
new_cols<-c("Year", stoplight$`Ecosystem Indicators`)
#switch rows and columns
stoplight2<-as_tibble(cbind(nms = names(stoplight), t(stoplight))) %>% filter(nms != "Ecosystem Indicators")
colnames(stoplight2)<-new_cols
stoplight2<-as_tibble(sapply(stoplight2, as.numeric))
head(stoplight2)

#making temperature names easier:
temp_indicators<-c("SST", "Summer20m", "Winter20m", "DeepTemp")
colnames(stoplight2)[5:8]<-temp_indicators

#just temperature indicators
temps<-stoplight2 %>% dplyr::select("Year", "SST", "Summer20m", "Winter20m", "DeepTemp") %>%
  mutate_at(temp_indicators, ~c(scale(.x))) #standardize

temps$time<-temps$Year-min(temps$Year) + 1
colnames(temps)<-c("Year", temp_indicators, "time")
summary(temps)

#making long for plotting
temps_long<-temps %>% select(-time) %>% pivot_longer(-c(Year), names_to = "Indicator", values_to = "Value")
ggplot(temps_long) + geom_point(aes(x = Year, y = Value)) + facet_wrap(~Indicator, scales = "free")

#state color palette
pal<-c("#61599d", "#c36377", "#f2af4a")
# Fit HMM -----------------------------------------------------------------
#3-state model for the physical indicators

temp_dists<-lapply(temp_indicators, function(x){x = "norm"})
names(temp_dists)<-temp_indicators
temps_hid <- MarkovChain$new(data = temps, n_states = 3) #3 state model
#time effect on mean, different for every state
f1<-list(SST = list(mean = ~ time),
         Summer20m = list(mean = ~ time),
         Winter20m = list(mean = ~ time), 
         DeepTemp = list(mean = ~time))


#fixing the slope parameters so they are the same across states
fixpar<-list(obs = c("SST.mean.state1.time" = 1, 
                     "SST.mean.state2.time" = 1, 
                     "SST.mean.state3.time" = 1, 
                     "Summer20m.mean.state1.time" = 2, 
                     "Summer20m.mean.state2.time" = 2, 
                     "Summer20m.mean.state3.time" = 2, 
                     "Winter20m.mean.state1.time" = 3, 
                     "Winter20m.mean.state2.time" = 3,
                     "Winter20m.mean.state3.time" = 3,
                     "DeepTemp.mean.state1.time" = 4,
                     "DeepTemp.mean.state2.time" = 4,
                     "DeepTemp.mean.state3.time" = 4))


# We will reconstruct / refit the model a large number of times. Even with the same
# seed at the top of the loop, the results from fit() may differ slightly run to run,
# even with the same seeds, because (I think) of hmmTMB's call to optimx() for optimization
AICs <- 0
best_aic <- 1.0e10

for(i in 1:500){
  print(i)
  set.seed(i)
  temps_inits2 <- temps_obs2 <- temps_hmm2 <- out1 <- NULL;
  # inits are generated here to bracket the extremes / labeling, and also include random component
  temps_inits2<-list()
  
  inits <- kmeans(temps[,2:5], centers = 3)$centers
  
  temps_inits2$SST<-list(mean = inits[,"SST"] + runif(3,-0.05,0.05), sd= rep(runif(1,0.8,1.2),3))
  temps_inits2$Summer20m<-list(mean = inits[,"Summer20m"] + runif(3,-0.05,0.05), sd= rep(runif(1,0.8,1.2),3))
  temps_inits2$Winter20m<-list(mean = inits[,"Winter20m"] + runif(3,-0.05,0.05), sd= rep(runif(1,0.8,1.2),3))
  temps_inits2$DeepTemp<-list(mean = inits[,"DeepTemp"] + runif(3,-0.05,0.05), sd= rep(runif(1,0.8,1.2),3))
  
  
  temps_obs2 <- Observation$new(data = temps, n_states = 3, dists = temp_dists, par = temps_inits2, formula = f1)
  temps_hmm2<-HMM$new(obs = temps_obs2, hid = temps_hid, fixpar = fixpar)
  
  for(t in 1:10){ #sometimes running the model multiple times with same initial values gets it to converge
    temps_hmm2$fit(silent = TRUE, itnmax = 2000)
    out1<-temps_hmm2$out() #save optimizer output to get convergence
    if(out1$convcode < 1){
      break
    }
  }
  
  if(out1$convcode < 1){ #checking convergence, 0 indicates successful 
    AICs[i] <- temps_hmm2$AIC_conditional()

    if(AICs[i] < best_aic) {
      best_model <- temps_hmm2
      best_aic <- AICs[i]
      print(best_aic)
      print(best_aic)
      #saved_list <- list(inits = temps_inits2, temps_obs2 = temps_obs2, temps_hmm2 = temps_hmm2, best_model = emps_hmm2)
    }
  }
}

temps_hmm2 <- best_model # use the best model across runs
temps_hmm2$out() 

#saving model object for later
#save(temps_hmm2, file = "Code/Results/physical_indicators_hmm.RData")

temps_hmm2$viterbi()
temps_hmm2$coeff_fe()
temps_hmm2$AIC_conditional()




# Manuscript Plots --------------------------------------------------------
#read in model fits if needed: 
load("Results/physical_indicators_hmm.RData")

#function to make dataframes out of parameter results=
Summ_table<-function(obspar, nstate ){
  obspar1<-obspar %>% as_tibble(rownames = "Ind") %>% 
    separate_wider_delim(Ind,delim = ".", names = c("Indicator", "par")) %>% 
    pivot_longer(all_of(3:(2+nstate)), names_to = "state", values_to = "estimate") %>%
    pivot_wider(id_cols = c(Indicator, state), names_from = "par", values_from = "estimate") %>% 
    mutate(state = as.integer(str_remove_all(state, "state ")))
  return(obspar1)
}

#making State 1 low, state 3 high, and state 2 mid
#right now state 1 is high, state 2 is low and state 3 is mid
#State 1 becomes state 3
#State 2 becomes state 1
#state 3 becomes state 2
#plotting like seabirds--model estimated means and confidence intervals plus data
#getting real values of temperatures
temp_means<-apply(stoplight2[, temp_indicators], 2, mean, na.rm = TRUE)
temp_sd<-apply(stoplight2[, temp_indicators], 2, sd, na.rm = TRUE)

#getting states
temps_sts<-tibble(Year = seq(1998, 2022, by = 1), state = temps_hmm2$viterbi()) %>% mutate(new_state = ifelse(state == 1, 3, ifelse(state == 2, 1, 2)))
temps_long_res<-temps_long %>% left_join(temps_sts, by = "Year") %>% 
  mutate(time = Year-min(Year) + 1) %>%
  mutate(Value_real = (Value * temp_sd[Indicator]) + temp_means[Indicator]) #adding in real values for plotting on real scale


#observation model parameters at each time step
obspar<-temps_hmm2$par(t = 1:25)
obs_ests<-apply(obspar$obs, 3, Summ_table, nstate = 3)
#combine
obs_ests2<-bind_rows(obs_ests, .id = "time") %>% mutate(time = as.numeric(time)) %>% 
  mutate(lower = qnorm(0.025, mean, sd), upper = qnorm(0.975, mean, sd)) %>%
  mutate(mean_real = (mean * temp_sd[Indicator]) + temp_means[Indicator], #adding in real values for plotting
         lower_real = (lower * temp_sd[Indicator]) + temp_means[Indicator],
         upper_real = (upper * temp_sd[Indicator]) + temp_means[Indicator])

#ordering states from low to high temperature (1 is low, 3 is high)
#state names stay the same
#color goes from low to high
names(pal)<-c("1", "3", "2")



temps_long_res<- temps_long_res %>%left_join(obs_ests2, by = c("time", "Indicator", "state"))

lab<-c("Sea surface\ntemperature (°C)", "Summer temperature\nat 20 m (°C)", "Winter temperature\nat 20 m (°C)", "Summer temperature\nat 50 m (°C)", "Summer salinity\nat 50 m (°C)")

names(lab)<-temp_indicators

ggplot(temps_long_res) + 
  geom_point(aes(x = Year, y = mean_real, color = as.factor(new_state), shape = as.factor(new_state))) +
  geom_linerange(aes(x = Year, ymin = lower_real, ymax = upper_real, color = as.factor(new_state))) +
  geom_point(aes(x = Year, y = Value_real)) + 
  scale_x_continuous(breaks = seq(1998, 2022, by = 4)) +
  facet_wrap(~factor(Indicator,levels = temp_indicators), strip.position = "left", 
             labeller = as_labeller(lab), nrow = 3) + 
  ylab(NULL) + 
  scale_color_manual(name = "State", values = pal) +
  scale_shape_manual(name = "State", values = c(8,17,15)) +
  theme_light() + theme(strip.background = element_blank(), strip.placement = "outside", 
                        strip.text = element_text(color = "black", size = 13), 
                        #legend.position = c(0.53,0.15), 
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 13),
                        legend.text = element_text(size = 12), 
                        legend.title = element_text(size =13), 
                        panel.grid.minor = element_blank())


ggsave("Figures/Physical_indicators_data_and_model_preds.png", dpi = 600)

#make time plot manually
#gives  parameter values at each time step
#from predict function (uncertainty around mean)
mod_preds<-temps_hmm2$predict("obspar", t=1:25, n_post = 1000)



#make summary table for each time step
estimates<-apply(mod_preds$mean, 3, Summ_table, nstate = 3)
dimnames(mod_preds$lcl)<-dimnames(mod_preds$mean)
dimnames(mod_preds$ucl)<-dimnames(mod_preds$mean)
estimateslwr<-apply(mod_preds$lcl, 3, Summ_table, nstate = 3)
estimatesupr<-apply(mod_preds$ucl, 3, Summ_table, nstate = 3)

#combine
mod_ests<-bind_rows(estimates, .id = "time") %>% mutate(time = as.numeric(time))
mod_lwr<-bind_rows(estimateslwr, .id = "time") %>% mutate(time = as.numeric(time)) %>% rename("lwr" = "mean")
mod_upr<-bind_rows(estimatesupr, .id = "time") %>% mutate(time = as.numeric(time)) %>% rename("upr" = "mean")

mod_ests<-mod_ests %>% left_join(mod_lwr, by = c("time", "Indicator", "state")) %>% 
  left_join(mod_upr, by = c("time", "Indicator", "state"))

#rename states
mod_ests1<-mod_ests %>% mutate(new_state = ifelse(state == 1, 3, ifelse(state == 2, 1, 2))) %>%
  mutate(Year = time + 1997)

#plot
plot2<-ggplot(mod_ests1) + 
  geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr, fill = as.factor(new_state), group = as.factor(new_state)), alpha = 0.4) +
  geom_line(aes(x = Year, y = mean, color = as.factor(new_state), group = as.factor(new_state))) + 
  facet_wrap(~factor(Indicator, levels = temp_indicators), scales = "free_y", labeller = as_labeller(lab), nrow = 3, strip.position = "left") + 
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  labs(y = NULL, x = "Year", fill = "State", color = "State") +
  scale_x_continuous(breaks = seq(1998, 2022, by = 4)) +
  theme_light() + 
  theme(strip.background = element_blank(), strip.placement = "outside", 
                        strip.text = element_text(color = "black", size = 13), 
                        #legend.position = c(0.53,0.15), 
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 13),
                        legend.text = element_text(size = 12), 
                        legend.title = element_text(size =13), 
                        panel.grid.minor = element_blank())
plot2

ggsave("Figures/time_effect.png", plot2, dpi = 600)  

temps_hmm2$coeff_fe()
#SST is 0.016
#summer 0.0012
#Winter 0. 040
#Deep is 0.056

temps_hmm2$par()$tpm
