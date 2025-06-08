library(boot)

output=readRDS("Results/model_output.Rdata")
#See a summary of the parameter estimates
output.sum <- output$BUGSoutput$summary%>%
  as.data.frame()%>%
  mutate(parameter=as.character(row.names(output.sum)))

#### Summarise occupancy

occupancy_data=output.sum%>%
  filter(str_detect(parameter, 'alpha'))%>%
  mutate(Sppcode=parse_number(parameter))%>%
  left_join(unique.data.frame(select(D,Species,Sppcode)),by="Sppcode")
occupancy_summary=  occupancy_data%>%
  mutate(mean=inv.logit(occupancy_data$mean))%>%
  mutate(sd=inv.logit(occupancy_data$sd))%>%
  mutate(X2.5.=inv.logit(occupancy_data$`2.5%`))%>%
  mutate(X50.=inv.logit(occupancy_data$`50%`))%>%
  mutate(X97.5.=inv.logit(occupancy_data$`97.5%`))%>%
  select(-`2.5%`:-`97.5%`)%>%
  relocate(Species)

#TO ESTIMATE GROUP-LEVEL HYPER-PARAMETERS:
#Define the occupancy covariate effects where mbeta is the community-level hyper-parameter, 
#gbeta is the group-level hyper-parameter, and sbeta is the species-specific parameter
mbeta <- output$BUGSoutput$sims.list$mbeta
gbeta <- output$BUGSoutput$sims.list$gbeta
sbeta <- output$BUGSoutput$sims.list$sbeta

#Define the covariates and the groups
covs <- colnames(X)
sizes <- c("Control","Carnivore","Herbivore")
#Create a data frame where the number of rows is equal to the number of covariates * the 
#number of groups
group <- data.frame(expand.grid(covs,sizes), matrix(NA,length(covs)*length(sizes),4))
colnames(group) <- c("Factor","Group","Mean","SD","LCI","UCI")
#Create a loop estimating the reference group values
for (a in 1:length(covs)){
  group[a,3:6] <- c(mean(mbeta[,a]),sd(mbeta[,a]),quantile(mbeta[,a],c(0.025,0.975)))
}
#Create a second loop estimating the other group values
for (a in 1:length(covs)){
  for (b in 1:(length(sizes)-1)){
    sims <- mbeta[,a] + gbeta[,((b-1)*dim(X)[2]+a)]
    group[(dim(X)[2]*(b)+a),3:6] <- c(mean(sims),sd(sims),quantile(sims,c(0.025,0.975)))
  }
}

#Export the results as a table
write.table(x=group,file="Results/group.csv",sep=",")

#TO ESTIMATE SPECIES-SPECIFIC COVARIATE VALUES:#Begin by defining the species
spec <- select(D,Species,Diet)%>%
  unique.data.frame()
#Define the group levels where 1 = Carnivore, 2 = herbivore, and 3 = omnivore
spec$Diet <- factor(spec$Diet)
gg <- as.numeric(spec$Diet)
#Define the occupancy covariates and groups
covs <- colnames(X)
sizes <- c("Omnivore","Carnivore","Herbivore")
#Create a data frame where the number of rows is equal to the number of covariates * the 
#number of species
species <- data.frame(expand.grid(covs,spec$Species), matrix(NA,length(covs)*length(spec$Species),dim(X)[2]))
colnames(species) <- c("Covariate","Species","Mean","SD","LCI","UCI")
#Re-define gbeta
gbeta <- output$BUGSoutput$sims.list$gbeta #Original
gbeta <- cbind(gbeta, matrix(0,nrow(gbeta),length(covs))) #New
#Create a loop that will estimate species-specific values for each of the covariates
for (a in 1:length(covs)){
  for (b in 1:length(spec$Species)){
    sims <- mbeta[,a] + gbeta[,((gg[b]-1)*dim(X)[2]+a)] + sbeta[,b,a]
    species[(dim(X)[2]*(b-1)+a),3:6] <- c(mean(sims), sd(sims),
                                          quantile(sims,c(0.025,0.975)))
  }
}

#Export the results as a table
write.csv(x=species,file="D:/WII-Thar_LTEO/Camera_trapping/MSOM/Trial/Results/species.csv",sep=",")

output=readRDS("D:/WII-Thar_LTEO/Camera_trapping/MSOM/Allmammals_output.Rdata")
#TO ESTIMATE SPECIES RICHNESS FOR EACH SITE:
#Begin by increasing memory to avoid errors
memory.limit(size=10000)
#Define the z matrix
z = output$BUGSoutput$sims.list$Z
#Sort the data frame based on species, study site, and diet category
d <- data.frame(ID = 1:nrow(D),D[,c("Species","Station","Diet")])#Create a new data frame
dz <- data.frame(d,t(z))
#Load reshape2 library
library(reshape2)
#Melt the data frame for easy casting
m.dz <- melt(dz,id.vars = c("Species","Station","Diet","ID") )
#Aggregate the data by summing the values in the z matrix for each camera station during each 
#iteration 
z.all <- acast(m.dz, Station ~ variable, fun.aggregate = sum) ### Add up 1/0 values for each species per iteration
#Use the aggregated values to create probability distributions and estimate mean, sd, and 95% 
#credible interval values for camera-station specific species richness
z.all <- t(apply(z.all,1,function(x) c(mean(x),sd(x),quantile(x,c(0.025,0.975))))); 
colnames(z.all) = c("Mean","SD","LCI","UCI")
#Export estimates of species richness as a table
write.table(x=z.all,file="Results/spprich.csv",sep=",")

#To estimate group richness for each site:
#Aggregate the data by summing the group-specific values in the z matrix for each camera 
#station during each iteration
z.group <- acast(m.dz,Station + Diet ~ variable, fun.aggregate = sum)
#Use the aggregated values to create probability distributions representing estimated camera-
#station specific group richness
z.group <- t(apply(z.group,1,function(x) c(mean(x),sd(x),quantile(x,c(0.025,0.975))))); 
colnames(z.group) = c("Mean","SD","LCI","UCI")
#Export estimates of group richness as a table
write.table(x=z.group,file="Results/grouprich.csv",sep=",")

left_join(data.frame(Sppcode=seq(1:16),psi.mean=inv.logit(output$BUGSoutput$mean$alpha),
                     psi.sd=inv.logit(output$BUGSoutput$sd$alpha)),
          data.frame(unique.data.frame(select(D,Species,Sppcode))),by="Sppcode")

naive_occu=pivot_wider(data = D,id_cols = Species,names_from = Station,values_from = detections,values_fn = sum)


