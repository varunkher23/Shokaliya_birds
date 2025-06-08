library(tidyverse)
library(readxl)

data = read.csv("Input/detection_history.csv")
survey_cov = read_xlsx("Input/survey_covariates.xlsx")%>%
  pivot_longer()
  
num.samples=data%>%group_by(Site)%>%
  summarise(num.samples=max(Survey))%>%ungroup()
Spp = data%>%select(Species)%>%unique.data.frame()%>%
  mutate(Sppcode=row_number())%>%
  cbind(Diet = rbinom(49,1,0.5))
sites = data%>%select(Site)%>%unique.data.frame()%>%
  cbind(sitecov1 = rbinom(34,1,0.5))

D = data%>%
  group_by(Species,Site)%>%
  summarise(detection = sum(Presence,na.rm = 1))%>%
  ungroup()%>%
  left_join(num.samples)%>%
  left_join(Spp)%>%
  left_join(sites)

#write.csv(D,"Input/MSOM_datafile.csv")

