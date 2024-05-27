#Code for running 
#Case Study 2 with Seabird density anomlaies 
#from Northern California current
library(tidyverse)
library(hmmTMB) #for HMM models
library(patchwork) #for plotting



# Data --------------------------------------------------------------------
#read data
NoCC_noFS<-read_csv("Data/NoCC_density_anomaly_noFS.csv")

#convert codes to species names
seabird_codes<-read_csv("Data/seabird_codes.csv")

#adding the NA years back in
a2<-seq(2003, 2022, 1)
NoCC_noFS<-NoCC_noFS %>% add_row(Year = a2[which(!(a2 %in% NoCC_noFS$Year))]) %>% arrange(Year)

#adding time column

NoCC_noFS<-NoCC_noFS %>% mutate(time = Year-min(Year) + 1)
summary(NoCC_noFS)


# Plotting Data -----------------------------------------------------------

#Easier to plot from long format
NoCC_noFS_long<-NoCC_noFS %>% pivot_longer(-c(Year,time), names_to = "Species", values_to = "Density_Anomaly")
ggplot(NoCC_noFS_long) + geom_point(aes(x = Year, y = Density_Anomaly)) + facet_wrap(~Species)

# Models ---------------------------------------------------------
#limiting species based on range of density anomaly values
NoCC_noFS_sd<-apply(NoCC_noFS, 2, sd,na.rm=T)[-c(1,56)]
NoCC_noFS_subset<-names(NoCC_noFS_sd)[which(NoCC_noFS_sd >= 0.2)]
#removing immature and unidentified gulls
idx<-which(NoCC_noFS_subset %in% c("IMGU", "GULL"))
NoCC_noFS_subset<-NoCC_noFS_subset[-idx]

#distributions
NoCC_dists<-lapply(NoCC_noFS_subset, function(x){x = "norm"})
names(NoCC_dists)<-sort(NoCC_noFS_subset)

#plotting species
NoCC_noFS_long %>% filter(Species %in% NoCC_noFS_subset) %>%
  ggplot() + geom_point(aes(x = Year, y = Density_Anomaly)) + facet_wrap(~Species, scales = "free_y")


#random initial values and model fitting
AICs <- 0
best_aic <- 1.0e10
best_model<-NULL
for(i in 1:200){
  print(i)
  NoCC_inits2_noFS<-NoCC_noFS_obs2<-NoCC_noFS_hmm2<-out1 <- NULL;
  
  #set up initial values
  NoCC_inits2_noFS<-list()
  inits<-kmeans(na.omit(NoCC_noFS[,NoCC_noFS_subset]), 2)$centers
  for(j in 1:length(NoCC_noFS_subset)){
    NoCC_inits2_noFS[[NoCC_noFS_subset[j]]]<-list(mean = inits[,j]+ runif(2,-0.05,0.05), sd= rep(runif(1,0.8,1.2),2))
  }
  
  #setting up model
  NoCC_noFS_hid2 <- MarkovChain$new(data = NoCC_noFS, n_states = 2)
  NoCC_noFS_obs2 <- Observation$new(data = NoCC_noFS, n_states = 2, dists = NoCC_dists, par = NoCC_inits2_noFS)
  
  NoCC_noFS_hmm2 <- HMM$new(obs = NoCC_noFS_obs2, hid = NoCC_noFS_hid2)
  
  #fit model
  NoCC_noFS_hmm2$fit(silent = TRUE)
  out<-NoCC_noFS_hmm2$out()
  
  if(out$convcode < 1){
    AICs[i] <- NoCC_noFS_hmm2$AIC_conditional()
    
    if(AICs[i] < best_aic) {
      best_model <- NoCC_noFS_hmm2
      best_aic <- AICs[i]
      print(best_aic)
    }
  }
}

NoCC_noFS_hmm2<-best_model
NoCC_noFS_hmm2$out()

NoCC_noFS_hmm2$par() 
NoCC_noFS_hmm2$viterbi() 


NoCC_noFS_hmm2$AIC_conditional()

#3-state
#random initial values and model fitting
AICs <- 0
best_aic <- 1.0e10
best_model<-NULL
for(i in 1:200){
  print(i)
  NoCC_inits3_noFS<-NoCC_noFS_obs3<-NoCC_noFS_hmm3<-out1 <- NULL;
  
  #set up initial values
  NoCC_inits3_noFS<-list()
  inits<-kmeans(na.omit(NoCC_noFS[,NoCC_noFS_subset]), 3)$centers
  for(j in 1:length(NoCC_noFS_subset)){
    NoCC_inits3_noFS[[NoCC_noFS_subset[j]]]<-list(mean = inits[,j]+ runif(3,-0.05,0.05), sd= rep(runif(1,0.8,1.2),3))
  }
  
  #setting up model
  NoCC_noFS_hid3 <- MarkovChain$new(data = NoCC_noFS, n_states = 3)
  NoCC_noFS_obs3 <- Observation$new(data = NoCC_noFS, n_states = 3, dists = NoCC_dists, par = NoCC_inits3_noFS)
  
  NoCC_noFS_hmm3 <- HMM$new(obs = NoCC_noFS_obs3, hid = NoCC_noFS_hid3)
  
  #fit model
  NoCC_noFS_hmm3$fit(silent = TRUE)
  out<-NoCC_noFS_hmm3$out()
  
  if(out$convcode < 1){
    AICs[i] <- NoCC_noFS_hmm3$AIC_conditional()
    
    if(AICs[i] < best_aic) {
      best_model <- NoCC_noFS_hmm3
      best_aic <- AICs[i]
      print(best_aic)
    }
  }
}


NoCC_noFS_hmm3<-best_model
NoCC_noFS_hmm3$par() 
NoCC_noFS_hmm3$viterbi() 


NoCC_noFS_hmm3$AIC_conditional()


# Results -----------------------------------------------------------------


#summary table
Summ_table<-function(obspar, nstate ){
  obspar1<-obspar %>% as_tibble(rownames = "Ind") %>% 
    separate_wider_delim(Ind,delim = ".", names = c("Indicator", "par")) %>% 
    pivot_longer(all_of(3:(2+nstate)), names_to = "state", values_to = "estimate") %>%
    pivot_wider(id_cols = c(Indicator, state), names_from = "par", values_from = "estimate") %>% 
    mutate(state = as.integer(str_remove_all(state, "state "))) 
  return(obspar1)
}

#NoCC
par2<-NoCC_noFS_hmm2$par()
obs_ests<-Summ_table(par2$obs, 2)
NoCC_noFS_tab<-obs_ests %>% 
  mutate(lower = qnorm(0.025, mean, sd), upper = qnorm(0.975, mean, sd))


#write_csv(NoCC_noFS_tab, "Code/NoCC_noFS_results_tab.csv")

#estimated states
NoCC_noFS_sts<-tibble(Year = seq(2003, 2022, by = 1), state = NoCC_noFS_hmm2$viterbi())
NoCC_noFS_long_res<-NoCC_noFS_long %>% filter(Species %in% NoCC_noFS_subset) %>% left_join(NoCC_noFS_sts, by = "Year") %>% 
  left_join(NoCC_noFS_tab, by = c("Species" = "Indicator", "state")) %>% left_join(seabird_codes, by = c("Species" = "CODE"))

col<-c("#11c2b5", "#677e8e")
#plot data and estimated states
plot1<-NoCC_noFS_long_res %>%
  ggplot() + 
  geom_point(aes(x = Year, y = mean, color = as.factor(state), shape = as.factor(state))) +
  geom_linerange(aes(x = Year, ymin = lower, ymax = upper, color = as.factor(state))) +
  geom_point(aes(x = Year, y = Density_Anomaly)) + 
  scale_shape_manual(values = c(8,17)) +
  scale_color_manual(values = col) + 
  scale_x_continuous(breaks = seq(2003, 2022, by = 6)) +
  labs(y = "Density anomaly") + 
  facet_wrap(~`Common Name`, labeller = labeller(`Common Name` = label_wrap_gen(15))) + 
  theme_light() + theme(legend.position = "none", 
                        strip.background = element_blank(), 
                        strip.text = element_text(color = "black", size = 11), 
                        axis.text = element_text(size = 10),
                        axis.title = element_text(size = 13),
                        legend.text = element_text(size = 12), 
                        legend.title = element_text(size =13), 
                        panel.grid.minor = element_blank())

plot1
#ggsave("Figures/SeabirdDensityresults.png", plot1, dpi = 600)

