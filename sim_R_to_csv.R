
load("W:/DATA/markscan_authorized_users/bergquist/EI/df_sim_full.RData")
load("W:/DATA/markscan_authorized_users/bergquist/EI/df_sim_250.RData")


write.csv(df_sim_250, file="C:/Users/bergquis/Dropbox/Evil_Insurers_Draft/unprofits/df_sim_250.csv")
write.csv(df_sim, file="C:/Users/bergquis/Dropbox/Evil_Insurers_Draft/unprofits/df_sim.csv")
