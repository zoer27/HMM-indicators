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
NoCC_noFS_subset<-names(NoCC_noFS_sd)[which(NoCC_noFS_sd >= 0.2)]
#removing immature and unidentified gulls
idx<-which(NoCC_noFS_subset %in% c("IMGU", "GULL"))
NoCC_noFS_subset<-NoCC_noFS_subset[-idx]

#distributions
NoCC_dists<-lapply(NoCC_noFS_subset, function(x){x = "norm"})
names(NoCC_dists)<-sort(NoCC_noFS_subset)

#initial values
#based on plots and what makes it converge
NoCC_noFS_long %>% filter(Species %in% NoCC_noFS_subset) %>%
  ggplot() + geom_point(aes(x = Year, y = Density_Anomaly)) + facet_wrap(~Species, scales = "free_y")

NoCC_inits2_noFS<-list()

NoCC_inits2_noFS$BFAL<-list(mean = c(0.5, -0.1), sd = c(1,1))
NoCC_inits2_noFS$CAAU<-list(mean = c(-0.3, 0.5), sd = c(1,1))
NoCC_inits2_noFS$CATE<-list(mean = c(0.25, -0.2), sd = c(1,1))
NoCC_inits2_noFS$COMU<-list(mean = c(0, 0.3), sd = c(1,1))
#NoCC_inits2_noFS$GULL<-list(mean = c(0.2, -0.1), sd = c(1,1))
#NoCC_inits2_noFS$IMGU<-list(mean = c(0.3, -0.1), sd = c(1,1))
NoCC_inits2_noFS$NOFU<-list(mean = c(0.4, -0.1), sd = c(1,1))
NoCC_inits2_noFS$PFSH<-list(mean = c(-0.1, 0.3), sd = c(1,1))
NoCC_inits2_noFS$SOSH<-list(mean = c(0.3, -0.8), sd = c(1,1))
NoCC_inits2_noFS$WEGU<-list(mean = c(0.4, -0.1), sd = c(1,1))
NoCC_inits2_noFS$WGGU<-list(mean = c(0.3, -0.1), sd = c(1,1))


#setting up model
NoCC_noFS_hid2 <- MarkovChain$new(data = NoCC_noFS, n_states = 2)
NoCC_noFS_obs2 <- Observation$new(data = NoCC_noFS, n_states = 2, dists = NoCC_dists, par = NoCC_inits2_noFS)

NoCC_noFS_hmm2 <- HMM$new(obs = NoCC_noFS_obs2, hid = NoCC_noFS_hid2)

#fit model
NoCC_noFS_hmm2$fit(silent = TRUE)
NoCC_noFS_hmm2$out()

NoCC_noFS_hmm2$par() 
NoCC_noFS_hmm2$viterbi() 

NoCC_noFS_hmm2$AIC_conditional()

#3-state
NoCC_inits3_noFS<-list()

NoCC_inits3_noFS$BFAL<-list(mean = c(0.5,-0.3, -0.1), sd = c(1,1,1))
NoCC_inits3_noFS$CAAU<-list(mean = c(-0.3, 0.5, 0.2), sd = c(1,1,1))
NoCC_inits3_noFS$CATE<-list(mean = c(0.5, 0.1, -0.2), sd = c(1,1,1))
NoCC_inits3_noFS$COMU<-list(mean = c(0, 0.8, 0.3), sd = c(1,1,1))
#NoCC_inits3_noFS$GULL<-list(mean = c(0.2, 0.4, -0.1), sd = c(1,1,1))
#NoCC_inits3_noFS$IMGU<-list(mean = c(0.3, -0.25, -0.1), sd = c(1,1,1))
NoCC_inits3_noFS$NOFU<-list(mean = c(0.4, -0.1, 0), sd = c(1,1,1))
NoCC_inits3_noFS$PFSH<-list(mean = c(-0.2, 0.3, -0.1), sd = c(1,1,1))
NoCC_inits3_noFS$SOSH<-list(mean = c(0.3, -0.8, 0.5), sd = c(1,1,1))
NoCC_inits3_noFS$WEGU<-list(mean = c(0.4,0, -0.1), sd = c(1,1,1))
NoCC_inits3_noFS$WGGU<-list(mean = c(-0.5, 0.3, -0.1), sd = c(1,1,1))

#setting up model
NoCC_noFS_hid3 <- MarkovChain$new(data = NoCC_noFS, n_states = 3)
NoCC_noFS_obs3 <- Observation$new(data = NoCC_noFS, n_states = 3, dists = NoCC_dists, par = NoCC_inits3_noFS)

NoCC_noFS_hmm3 <- HMM$new(obs = NoCC_noFS_obs3, hid = NoCC_noFS_hid3)

#fit model
NoCC_noFS_hmm3$fit(silent = TRUE)
NoCC_noFS_hmm3$out()

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


#write_csv(NoCC_noFS_tab, "Code/Results/NoCC_noFS_results_tab.csv")

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

