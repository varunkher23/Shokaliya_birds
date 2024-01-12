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
Spp = data%>%select(Species)%>%unique.data.frame()%>%
  mutate(Sppcode=row_number())

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
                  det.formula = ~scale(day) + scale(tod), 
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
hist(prediction_occu$psi.0.samples[,17,1]) ### Scrubland
hist(prediction_occu$psi.0.samples[,17,34]) ### Cropland

prediction_det = predict(out.ms,cbind(1,survey_cov%>%select(Date,Time)%>%unique.data.frame()%>%scale()),type = "detection")
hist(prediction_det$p.0.samples[,12,10])