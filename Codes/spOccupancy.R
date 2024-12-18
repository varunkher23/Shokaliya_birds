#### spOccupancy

library(tidyverse)
library(readxl)
library(spOccupancy)
library(boot)

data = read.csv("Input/detection_history.csv")
survey_cov = read_xlsx("Input/survey_covariates.xlsx")%>%
  pivot_longer(`Start Time` : `10.0`,names_to = "Survey",values_to =  "End_Time")%>%
  filter(Survey != "Start Time")%>%
  mutate(Survey = as.numeric(Survey))%>%
  mutate(Date = yday(as.Date(End_Time)))%>%
  mutate(Time = as.numeric(difftime(End_Time, as.Date(End_Time),units = "mins"))-5)%>%
  arrange(Date)
sites = read_xlsx("Input/site_covariates.xlsx")  


num.samples=data%>%group_by(Site)%>%
  summarise(num.samples=max(Survey))%>%ungroup()
sp_cov = read_xlsx("Input/species_covariates.xlsx")%>%
  rename(Species=`Latin Name`)%>%
  rename(feeding = `Feeding Guild`)
Spp = data%>%select(Species)%>%unique.data.frame()%>%
  mutate(Sppcode=row_number())%>%
  left_join(sp_cov)

Spp%>%filter(feeding=="Granivorous")%>%select(Sppcode)%>%unlist%>%as.vector()

y = array(NA,dim = c(nrow(Spp),nrow(sites),max(num.samples$num.samples)))
dimnames(y)[[1]] <- Spp$Species
dimnames(y)[[2]] <- sites$Site
str(y)

for (j in 1:nrow(sites)){
  for (k in 1:max(num.samples$num.samples)) {
    curr.df = data%>%
      filter(Site == sites$Site[j], Survey == k)
    
    curr.sp = which(Spp$Species %in% (curr.df%>%filter(Presence == 1))$Species)
    y[curr.sp,j,k] = 1
    y[-curr.sp,j,k] = 0
  }
}

### Total observations for each species
# Divide by num.sites x num.surveys to get naive occupancy
apply(y, 1, sum, na.rm = TRUE)


####### Survey Covariates #######
hb.day <- matrix(NA, nrow = nrow(sites), max(num.samples$num.samples))
hb.tod <- matrix(NA, nrow = nrow(sites), max(num.samples$num.samples))

for (j in 1:nrow(sites)){
  for (k in 1:max(num.samples$num.samples)) {
    current.vals = survey_cov%>%
      filter(Site == sites$Site[j], Survey == k)
    hb.day[j,k] = current.vals$Date[1]
    hb.tod[j,k] = current.vals$Time[1]
  }
}

habitat = as.numeric(sites$`Land Cover` == "Scrub")
occ.covs = data.frame(habitat)

det.covs = list(day= hb.day, tod = hb.tod, habitat = habitat)

data.msom = list(y = y, occ.covs = occ.covs,
                 det.covs= det.covs)

#### Intitials

N = nrow(Spp)
ms.inits <- list(alpha.comm = 0,
                 beta.comm = 0,
                 beta = 0,
                 alpha = 0,
                 tau.sq.beta = 1,
                 tau.sq.alpha = 1,
                 z = apply(y,c(1,2),max,na.rm = T))

ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

### Run the model

# Approx. run time:  6 min
out.ms <- msPGOcc(occ.formula = ~habitat, 
                  det.formula = ~scale(day) + scale(tod) + habitat, 
                  data = data.msom, 
                  inits = ms.inits, 
                  n.samples = 30000, 
                  priors = ms.priors, 
                  n.omp.threads = 6, 
                  verbose = TRUE, 
                  n.report = 6000, 
                  n.burn = 10000,
                  n.thin = 50, 
                  n.chains = 3)

summary(out.ms, level = "community")

#### Predict
prediction_occu = predict(out.ms,cbind(1,habitat))

prediction_det = predict(out.ms,cbind(1,survey_cov%>%select(Date,Time)%>%unique.data.frame()%>%scale()),type = "detection")
prediction_det = predict(out.ms,cbind(1,mean(scale(survey_cov$Date)),mean(scale(survey_cov$Time)),c(1,0)),type = "detection")

species = "Prinia socialis"
hist(prediction_occu$psi.0.samples[,which(Spp$Species==species),1],xlim = c(0,1),col = "yellow") ### Scrubland
hist(prediction_occu$psi.0.samples[,which(Spp$Species==species),34],xlim = c(0,1),col = "green",add=T) ### Cropland
hist(prediction_det$p.0.samples[,which(Spp$Species==species),1],xlim = c(0,1),col = "yellow")
hist(prediction_det$p.0.samples[,which(Spp$Species==species),2],xlim = c(0,1),col = "green",add=T)

### Richness
rich.samples <- apply(prediction_occu$z.0.samples, c(1, 3), sum)
rich.median <- apply(rich.samples, 2, median, na.rm = TRUE)
rich.low <- apply(rich.samples, 2, quantile, 0.025, na.rm = TRUE)
rich.high <- apply(rich.samples, 2, quantile, 0.975, na.rm = TRUE)
richness.ci.width = rich.high - rich.low
richness = data.frame(median = rich.median,lcl = rich.low, ucl = rich.high, habitat = habitat)%>%unique.data.frame()




#### Function to get habitat specific occupancy for species
get_sp_occupancy <- function(species){
  curr.sp <- which(Spp$Species == species)
  curr.sp.psi.samples <- prediction_occu$psi.0.samples[, curr.sp, ]
  curr.sp.occ <- apply(curr.sp.psi.samples, 2, mean)
  curr.sp.occ.lcl <- apply(curr.sp.psi.samples, 2, quantile, 0.025)
  curr.sp.occ.ucl <- apply(curr.sp.psi.samples, 2, quantile, 0.975)
  curr.sp.occ.ci.width = curr.sp.occ.ucl - curr.sp.occ.lcl
  sp.occ = data.frame(Species = species,occ_median = curr.sp.occ, occ_lcl = curr.sp.occ.lcl, occ_ucl = curr.sp.occ.ucl, curr.sp.ci.width = curr.sp.occ.ci.width,
                      habitat = habitat)%>%unique.data.frame()%>%mutate(habitat = ifelse(habitat==as.numeric(1),"Scrub","Crop"))
  output = list(sp.occ = sp.occ, sp.psi.samples = curr.sp.psi.samples)
  return(output)
}

sp_occupancy_temp = get_sp_occupancy(species = 'Coturnix coromandelica')$sp.occ

### Guildwise species richness

posterior_guildwise = data.frame() 
for (i in 1:length(unique(sp_cov$feeding))) {
guild = Spp%>%filter(feeding==sp_cov$feeding[i])%>%select(Sppcode)%>%unlist%>%as.vector()
rich.guild = apply(prediction_occu$z.0.samples[,guild,], c(1, 3), sum)
posterior_guild = rbind(data.frame(sprich = rich.guild[,1],Habitat = "Scrubland",feeding = sp_cov$feeding[i]),
                            data.frame(sprich = rich.guild[,17],Habitat = "Cropland", feeding = sp_cov$feeding[i]))
posterior_guildwise = rbind(posterior_guildwise,posterior_guild)
}


library(ggbeeswarm)
ggplot(posterior_guildwise%>%filter(feeding%in%c("Insectivorous","Granivorous", "Omnivorous")), 
       aes(x = Habitat, y = sprich, col = Habitat)) +
  geom_quasirandom()+
  theme_classic()+
  ylim(c(0,15))+
  facet_grid(~feeding)

