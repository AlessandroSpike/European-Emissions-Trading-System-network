library(memisc)
library(reshape2)
library(stats)
library(Hmisc)
library(texreg)
library(xtable)
library(ggplot2)
library(plm)
library(ggthemes)
library(grid)
library(gridExtra)
library("scales")
library(gtable)
library(plyr)
library(statnet)
library(igraph)
library(rvest)
library(readr)
library(lubridate)
library(data.table)
library(stringdist)
library(dplyr)
library(plm)
library(lmtest)
library(Hmisc)
library(stargazer)
library(grid)
library(gridExtra)
library(egg)
library(pglm)
library(stats)
library(tsibble)
library(giscoR)
library(gravity)
library("RColorBrewer")
library(viridis)
library(colorspace)
library(scico)



#############
# IMPORT DB #
#############

db_model <- read_csv("../db_merge_final3.csv")[,-1]

# variable same country 
db_model$country_iso_o <- substring(db_model$iso_o,1,2)
db_model$country_iso_d <- substring(db_model$iso_d,1,2)
db_model$country_same <- ifelse(db_model$country_iso_o==db_model$country_iso_d,1,0)

# variable same sector 
db_model$same_sector_2005 <- ifelse(db_model$sector_d.2005==db_model$sector_o.2005,1,0)
db_model$same_sector_2006 <- ifelse(db_model$sector_d.2006==db_model$sector_o.2006,1,0)
db_model$same_sector_2007 <- ifelse(db_model$sector_d.2007==db_model$sector_o.2007,1,0)
db_model$same_sector_2008 <- ifelse(db_model$sector_d.2008==db_model$sector_o.2008,1,0)
db_model$same_sector_2009 <- ifelse(db_model$sector_d.2009==db_model$sector_o.2009,1,0)
db_model$same_sector_2010 <- ifelse(db_model$sector_d.2010==db_model$sector_o.2010,1,0)
db_model$same_sector_2011 <- ifelse(db_model$sector_d.2011==db_model$sector_o.2011,1,0)
db_model$same_sector_2012 <- ifelse(db_model$sector_d.2012==db_model$sector_o.2012,1,0)
db_model$same_sector_2013 <- ifelse(db_model$sector_d.2013==db_model$sector_o.2013,1,0)
db_model$same_sector_2014 <- ifelse(db_model$sector_d.2014==db_model$sector_o.2014,1,0)
db_model$same_sector_2015 <- ifelse(db_model$sector_d.2015==db_model$sector_o.2015,1,0)
db_model$same_sector_2016 <- ifelse(db_model$sector_d.2016==db_model$sector_o.2016,1,0)
db_model$same_sector_2017 <- ifelse(db_model$sector_d.2017==db_model$sector_o.2017,1,0)
db_model$same_sector_2018 <- ifelse(db_model$sector_d.2018==db_model$sector_o.2018,1,0)
db_model$same_sector_2019 <- ifelse(db_model$sector_d.2019==db_model$sector_o.2019,1,0)
db_model$same_sector_2020 <- ifelse(db_model$sector_d.2020==db_model$sector_o.2020,1,0)


# YEAR: 2005
fit_2005 <- ppml(dependent_variable = "amount.2005", distance = "distance",
                additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                "POP_2005_acq", "POP_2005_transf", 
                "GDP_2005_acq", "GDP_2005_transf",
                "allocation.2005_acq","allocation.2005_transf",
                "emissions.2005_acq","emissions.2005_transf","country_same","same_sector_2005"),
                robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2006
fit_2006 <- ppml(dependent_variable = "amount.2006", distance = "distance",
                additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                "POP_2006_acq", "POP_2006_transf", 
                "GDP_2006_acq", "GDP_2006_transf",
                "allocation.2006_acq","allocation.2006_transf",
                "emissions.2006_acq","emissions.2006_transf","country_same","same_sector_2006"),
                robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2007
fit_2007 <- ppml(dependent_variable = "amount.2007", distance = "distance",
                additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                "POP_2007_acq", "POP_2007_transf", 
                "GDP_2007_acq", "GDP_2007_transf",
                "allocation.2007_acq","allocation.2007_transf",
                "emissions.2007_acq","emissions.2007_transf","country_same","same_sector_2007"),
                robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2008
fit_2008 <- ppml(dependent_variable = "amount.2008", distance = "distance",
                additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                "POP_2008_acq", "POP_2008_transf", 
                "GDP_2008_acq", "GDP_2008_transf",
                "allocation.2008_acq","allocation.2008_transf",
                "emissions.2008_acq","emissions.2008_transf","country_same","same_sector_2008"),
                robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2009
fit_2009 <- ppml(dependent_variable = "amount.2009", distance = "distance",
                additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                "POP_2009_acq", "POP_2009_transf", 
                "GDP_2009_acq", "GDP_2009_transf",
                "allocation.2009_acq","allocation.2009_transf",
                "emissions.2009_acq","emissions.2009_transf","country_same","same_sector_2009"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2010
fit_2010 <- ppml(dependent_variable = "amount.2010", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2010_acq", "POP_2010_transf", 
                 "GDP_2010_acq", "GDP_2010_transf",
                 "allocation.2010_acq","allocation.2010_transf",
                 "emissions.2010_acq","emissions.2010_transf","country_same","same_sector_2010"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2011
fit_2011 <- ppml(dependent_variable = "amount.2011", distance = "distance",
                additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                "POP_2011_acq", "POP_2011_transf", 
                "GDP_2011_acq", "GDP_2011_transf",
                "allocation.2011_acq","allocation.2011_transf",
                "emissions.2011_acq","emissions.2011_transf","country_same","same_sector_2011"),
                robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2012
fit_2012 <- ppml(dependent_variable = "amount.2012", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2012_acq", "POP_2012_transf", 
                 "GDP_2012_acq", "GDP_2012_transf",
                 "allocation.2012_acq","allocation.2012_transf",
                 "emissions.2012_acq","emissions.2012_transf","country_same","same_sector_2012"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2013
fit_2013 <- ppml(dependent_variable = "amount.2013", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2013_acq", "POP_2013_transf", 
                 "GDP_2013_acq", "GDP_2013_transf",
                 "allocation.2013_acq","allocation.2013_transf",
                 "emissions.2013_acq","emissions.2013_transf","country_same","same_sector_2013"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2014
fit_2014 <- ppml(dependent_variable = "amount.2014", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2014_acq", "POP_2014_transf", 
                 "GDP_2014_acq", "GDP_2014_transf",
                "allocation.2014_acq","allocation.2014_transf",
                 "emissions.2014_acq","emissions.2014_transf","country_same","same_sector_2014"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2015
fit_2015 <- ppml(dependent_variable = "amount.2015", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2015_acq", "POP_2015_transf", 
                 "GDP_2015_acq", "GDP_2015_transf",
                 "allocation.2015_acq","allocation.2015_transf",
                 "emissions.2015_acq","emissions.2015_transf","country_same","same_sector_2015"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2016
fit_2016 <- ppml(dependent_variable = "amount.2016", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2016_acq", "POP_2016_transf", 
                 "GDP_2016_acq", "GDP_2016_transf",
                 "allocation.2016_acq","allocation.2016_transf",
                 "emissions.2016_acq","emissions.2016_transf","country_same","same_sector_2016"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2017
fit_2017 <- ppml(dependent_variable = "amount.2017", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2017_acq", "POP_2017_transf", 
                 "GDP_2017_acq", "GDP_2017_transf",
                 "allocation.2017_acq","allocation.2017_transf",
                 "emissions.2017_acq","emissions.2017_transf","country_same","same_sector_2017"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2018
fit_2018 <- ppml(dependent_variable = "amount.2018", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2018_acq", "POP_2018_transf", 
                 "GDP_2018_acq", "GDP_2018_transf",
                 "allocation.2018_acq","allocation.2018_transf",
                 "emissions.2018_acq","emissions.2018_transf","country_same","same_sector_2018"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2019
fit_2019 <- ppml(dependent_variable = "amount.2019", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2019_acq", "POP_2019_transf", 
                 "GDP_2019_acq", "GDP_2019_transf",
                 "allocation.2019_acq","allocation.2019_transf",
                 "emissions.2019_acq","emissions.2019_transf","country_same","same_sector_2019"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")

# YEAR: 2020
fit_2020 <- ppml(dependent_variable = "amount.2020", distance = "distance",
                 additional_regressors = c("Total_area_SqKm_acq", "Total_area_SqKm_transf",
                 "POP_2020_acq", "POP_2020_transf", 
                 "GDP_2020_acq", "GDP_2020_transf",
                 "allocation.2020_acq","allocation.2020_transf",
                 "emissions.2020_acq","emissions.2020_transf","country_same","same_sector_2020"),
                 robust = T, method = "white1", data = db_model, cluster = "distance")


# final Table
stats_db <- data.frame()
stats_db_plot <- data.frame()
stats_db_plot_SE <- data.frame()
list_variables <- c("Intercept","Distance","Total_area_SqKm_acq", "Total_area_SqKm_transf",
                    "pop_acq", "pop_transf", 
                    "gdp_acq", "gdp_transf",
                    "allocation_acq","allocation_transf",
                    "emissions_acq","emissions_transf","country_same")


for (i in c(1:length(list_variables))){
  stats_db[i,1] <- list_variables[i]
}
names(stats_db)[1] <- "Variable" 

for (i in c(1:length(list_variables))){
  stats_db[i,2] <- paste0(round(fit_2005$coefficients[i,1], digits=3)," (",round(fit_2005$coefficients[i,2], digits=3),")")
}
names(stats_db)[2] <- "2005" 

for (i in c(1:length(list_variables))){
  stats_db[i,3] <- paste0(round(fit_2006$coefficients[i,1], digits=3)," (",round(fit_2006$coefficients[i,2], digits=3),")")
}
names(stats_db)[3] <- "2006" 

for (i in c(1:length(list_variables))){
  stats_db[i,4] <- paste0(round(fit_2007$coefficients[i,1], digits=3)," (",round(fit_2007$coefficients[i,2], digits=3),")")
}
names(stats_db)[4] <- "2007" 

for (i in c(1:length(list_variables))){
  stats_db[i,5] <- paste0(round(fit_2008$coefficients[i,1], digits=3)," (",round(fit_2008$coefficients[i,2], digits=3),")")
}
names(stats_db)[5] <- "2008" 

for (i in c(1:length(list_variables))){
  stats_db[i,6] <- paste0(round(fit_2009$coefficients[i,1], digits=3)," (",round(fit_2009$coefficients[i,2], digits=3),")")
}
names(stats_db)[6] <- "2009" 

for (i in c(1:length(list_variables))){
  stats_db[i,7] <- paste0(round(fit_2010$coefficients[i,1], digits=3)," (",round(fit_2010$coefficients[i,2], digits=3),")")
}
names(stats_db)[7] <- "2010" 

for (i in c(1:length(list_variables))){
  stats_db[i,8] <- paste0(round(fit_2011$coefficients[i,1], digits=3)," (",round(fit_2011$coefficients[i,2], digits=3),")")
}
names(stats_db)[8] <- "2011" 

for (i in c(1:length(list_variables))){
  stats_db[i,9] <- paste0(round(fit_2012$coefficients[i,1], digits=3)," (",round(fit_2012$coefficients[i,2], digits=3),")")
}
names(stats_db)[9] <- "2012" 

for (i in c(1:length(list_variables))){
  stats_db[i,10] <- paste0(round(fit_2013$coefficients[i,1], digits=3)," (",round(fit_2013$coefficients[i,2], digits=3),")")
}
names(stats_db)[10] <- "2013" 

for (i in c(1:length(list_variables))){
  stats_db[i,11] <- paste0(round(fit_2014$coefficients[i,1], digits=3)," (",round(fit_2014$coefficients[i,2], digits=3),")")
}
names(stats_db)[11] <- "2014" 

for (i in c(1:length(list_variables))){
  stats_db[i,12] <- paste0(round(fit_2015$coefficients[i,1], digits=3)," (",round(fit_2015$coefficients[i,2], digits=3),")")
}
names(stats_db)[12] <- "2015" 

for (i in c(1:length(list_variables))){
  stats_db[i,13] <- paste0(round(fit_2016$coefficients[i,1], digits=3)," (",round(fit_2016$coefficients[i,2], digits=3),")")
}
names(stats_db)[13] <- "2016" 

for (i in c(1:length(list_variables))){
  stats_db[i,14] <- paste0(round(fit_2017$coefficients[i,1], digits=3)," (",round(fit_2017$coefficients[i,2], digits=3),")")
}
names(stats_db)[14] <- "2017" 

for (i in c(1:length(list_variables))){
  stats_db[i,15] <- paste0(round(fit_2018$coefficients[i,1], digits=3)," (",round(fit_2018$coefficients[i,2], digits=3),")")
}
names(stats_db)[15] <- "2018" 

for (i in c(1:length(list_variables))){
  stats_db[i,16] <- paste0(round(fit_2019$coefficients[i,1], digits=3)," (",round(fit_2019$coefficients[i,2], digits=3),")")
}
names(stats_db)[16] <- "2019" 

for (i in c(1:length(list_variables))){
  stats_db[i,17] <- paste0(round(fit_2020$coefficients[i,1], digits=3)," (",round(fit_2020$coefficients[i,2], digits=3),")")
}
names(stats_db)[17] <- "2020" 

stats_db[length(list_variables)+1,1] <- "sector_same"
stats_db[length(list_variables)+1,2] <- paste0(round(fit_2005$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2005$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,3] <- paste0(round(fit_2006$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2006$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,4] <- paste0(round(fit_2007$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2007$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,5] <- paste0(round(fit_2008$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2008$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,6] <- paste0(round(fit_2009$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2009$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,7] <- paste0(round(fit_2010$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2010$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,8] <- paste0(round(fit_2011$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2011$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,9] <- paste0(round(fit_2012$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2012$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,10] <- paste0(round(fit_2013$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2013$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,11] <- paste0(round(fit_2014$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2014$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,12] <- paste0(round(fit_2015$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2015$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,13] <- paste0(round(fit_2016$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2016$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,14] <- paste0(round(fit_2017$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2017$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,15] <- paste0(round(fit_2018$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2018$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,16] <- paste0(round(fit_2019$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2019$coefficients[length(list_variables)+1,2], digits=3),")")
stats_db[length(list_variables)+1,17] <- paste0(round(fit_2020$coefficients[length(list_variables)+1,1], digits=3)," (",round(fit_2020$coefficients[length(list_variables)+1,2], digits=3),")")


# table no SE

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,1] <- list_variables[i]
}
names(stats_db_plot)[1] <- "Variable" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,2] <- round(fit_2005$coefficients[i,1], digits=20)
}
names(stats_db_plot)[2] <- "2005" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,3] <- round(fit_2006$coefficients[i,1], digits=20)
}
names(stats_db_plot)[3] <- "2006" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,4] <- round(fit_2007$coefficients[i,1], digits=20)
}
names(stats_db_plot)[4] <- "2007" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,5] <- round(fit_2008$coefficients[i,1], digits=20)
}
names(stats_db_plot)[5] <- "2008" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,6] <- round(fit_2009$coefficients[i,1], digits=20)
}
names(stats_db_plot)[6] <- "2009" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,7] <- round(fit_2010$coefficients[i,1], digits=20)
}
names(stats_db_plot)[7] <- "2010" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,8] <- round(fit_2011$coefficients[i,1], digits=20)
}
names(stats_db_plot)[8] <- "2011" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,9] <- round(fit_2012$coefficients[i,1], digits=20)
}
names(stats_db_plot)[9] <- "2012" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,10] <- round(fit_2013$coefficients[i,1], digits=20)
}
names(stats_db_plot)[10] <- "2013" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,11] <- round(fit_2014$coefficients[i,1], digits=20)
}
names(stats_db_plot)[11] <- "2014" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,12] <- round(fit_2015$coefficients[i,1], digits=20)
}
names(stats_db_plot)[12] <- "2015" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,13] <- round(fit_2016$coefficients[i,1], digits=20)
}
names(stats_db_plot)[13] <- "2016" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,14] <- round(fit_2017$coefficients[i,1], digits=20)
}
names(stats_db_plot)[14] <- "2017" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,15] <- round(fit_2018$coefficients[i,1], digits=20)
}
names(stats_db_plot)[15] <- "2018" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,16] <- round(fit_2019$coefficients[i,1], digits=20)
}
names(stats_db_plot)[16] <- "2019" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot[i,17] <- round(fit_2020$coefficients[i,1], digits=20)
}
names(stats_db_plot)[17] <- "2020" 

stats_db_plot[length(list_variables)+1,1] <- "sector_same"


# table ONLY SE

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,1] <- list_variables[i]
}
names(stats_db_plot_SE)[1] <- "Variable" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,2] <- round(fit_2005[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[2] <- "2005" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,3] <- round(fit_2006[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[3] <- "2006" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,4] <- round(fit_2007[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[4] <- "2007" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,5] <- round(fit_2008[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[5] <- "2008" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,6] <- round(fit_2009[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[6] <- "2009" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,7] <- round(fit_2010[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[7] <- "2010" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,8] <- round(fit_2011[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[8] <- "2011" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,9] <- round(fit_2012[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[9] <- "2012" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,10] <- round(fit_2013[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[10] <- "2013" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,11] <- round(fit_2014[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[11] <- "2014" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,12] <- round(fit_2015[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[12] <- "2015" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,13] <- round(fit_2016[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[13] <- "2016" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,14] <- round(fit_2017[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[14] <- "2017" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,15] <- round(fit_2018[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[15] <- "2018" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,16] <- round(fit_2019[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[16] <- "2019" 

for (i in c(1:nrow(stats_db))){
  stats_db_plot_SE[i,17] <- round(fit_2020[[1]][i,2], digits=20)
}
names(stats_db_plot_SE)[17] <- "2020" 

stats_db_plot_SE[length(list_variables)+1,1] <- "sector_same"



stats_db_plot_long <- reshape(data=stats_db_plot, 
        direction = "long",
        varying = 2:17,
        v.names = "Value",
        idvar = c("Variable"),
        timevar = "Year",
        times = c(2005:2020))
names(stats_db_plot_long)[3] <- "beta"

stats_db_plot_SE_long <- reshape(data=stats_db_plot_SE, 
                              direction = "long",
                              varying = 2:17,
                              v.names = "Value",
                              idvar = c("Variable"),
                              timevar = "Year",
                              times = c(2005:2020))
names(stats_db_plot_SE_long)[3] <- "S.E."



stats_db_plot_merge <- merge(stats_db_plot_long,stats_db_plot_SE_long, by=c("Variable","Year"))
stats_db_plot_merge$hi <- stats_db_plot_merge$beta + stats_db_plot_merge$S.E.
stats_db_plot_merge$lo <- stats_db_plot_merge$beta - stats_db_plot_merge$S.E.


plot_distance <- ggplot(aes(x=Year), 
                        data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="Distance"),])+
                        geom_line(aes(y=beta), color="blue",size=1)+
                        geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                        geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                        geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                        geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                        annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                        annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                        annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                        scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Distance") +xlab("")


plot_country <- ggplot(aes(x=Year), 
                        data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="country_same"),])+
                        geom_line(aes(y=beta), color="blue",size=1)+
                        geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                        geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                        geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                        geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                        annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                        annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                        annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                        scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: country_same") +xlab("")


plot_sector <- ggplot(aes(x=Year), 
                       data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="sector_same"),])+
                       geom_line(aes(y=beta), color="blue",size=1)+
                       geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                       geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                       geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                       geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                       annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                       annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                       annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                       scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: sector_same") +xlab("")

plot_allocation_acq <- ggplot(aes(x=Year),
                       data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="allocation_acq"),])+
                       geom_line(aes(y=beta), color="blue",size=1)+
                       geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                       geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                       geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                       geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                       annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                       annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                       annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                       scale_x_continuous(breaks=seq(2005,2020,5))+theme_bw() + ylab("Beta: Allocation purchaser")

plot_allocation_transf <- ggplot(aes(x=Year),
                          data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="allocation_transf"),])+
                          geom_line(aes(y=beta), color="blue",size=1)+
                          geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                          geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                          geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                          geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                          annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                          annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                          annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                          scale_x_continuous(breaks=seq(2005,2020,5))+theme_bw() + ylab("Beta: Allocation seller")


plot_emissions_acq <- ggplot(aes(x=Year), 
                      data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="emissions_acq"),])+
                      geom_line(aes(y=beta), color="blue",size=1)+
                      geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                      geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                      geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                      geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                      annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                      annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                      annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                      scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Emissions purchaser") +xlab("")

plot_emissions_transf <- ggplot(aes(x=Year), 
                         data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="emissions_transf"),])+
                         geom_line(aes(y=beta), color="blue",size=1)+
                         geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                         geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                         geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                         geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                         annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                         annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                         annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                         scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Emissions seller") +xlab("")

plot_pop_acq <- ggplot(aes(x=Year), 
                        data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="pop_acq"),])+
                        geom_line(aes(y=beta), color="blue",size=1)+
                        geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                        geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                        geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                        geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                        annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                        annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                        annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                        scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Population purchaser") +xlab("")

plot_pop_transf <- ggplot(aes(x=Year), 
                          data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="pop_transf"),])+
                          geom_line(aes(y=beta), color="blue",size=1)+
                          geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                          geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                          geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                          geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                          annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                          annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                          annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                          scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Population seller") +xlab("")


plot_area_acq <- ggplot(aes(x=Year), 
                        data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="Total_area_SqKm_acq"),])+
                        geom_line(aes(y=beta), color="blue",size=1)+
                        geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                        geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                        geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                        geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                        annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                        annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                        annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                        scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Total area SqKm purchaser") +xlab("")

plot_area_transf <- ggplot(aes(x=Year), 
                          data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="Total_area_SqKm_transf"),])+
                          geom_line(aes(y=beta), color="blue",size=1)+
                          geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                          geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                          geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                          geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                          annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                          annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                          annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                          scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: Total area SqKm seller") +xlab("")


plot_gdp_acq <- ggplot(aes(x=Year), 
                            data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="gdp_acq"),])+
                            geom_line(aes(y=beta), color="blue",size=1)+
                            geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                            geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                            geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                            geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                            annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                            annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                            annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                            scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: GDP purchaser") +xlab("")

plot_gdp_transf <- ggplot(aes(x=Year), 
                                data=stats_db_plot_merge[which(stats_db_plot_merge$Variable=="gdp_transf"),])+
                                geom_line(aes(y=beta), color="blue",size=1)+
                                geom_line(aes(y=hi), color="blue",size=0.5, linetype="dashed")+
                                geom_line(aes(y=lo), color="blue",size=0.5, linetype="dashed")+
                                geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
                                geom_vline(xintercept=c(2008,2013), linetype="dotted", color = "black", size=1)+
                                annotate('rect',xmin=2005,xmax=2008,ymin=-Inf,ymax=Inf,fill="green",alpha=0.05)+
                                annotate('rect',xmin=2008,xmax=2013,ymin=-Inf,ymax=Inf,fill="green",alpha=0.1)+
                                annotate('rect',xmin=2013,xmax=2020,ymin=-Inf,ymax=Inf,fill="green",alpha=0.15)+
                                scale_x_continuous(breaks=seq(2005,2020,2))+theme_bw() + ylab("Beta: GDP seller") +xlab("")
                              


plot_core <- ggarrange(plot_distance,
                      plot_country, plot_sector,
                      ncol = 1, nrow = 3)

plot_emission <- ggarrange(plot_emissions_acq, plot_emissions_transf, ncol = 1, nrow = 2)
plot_allocation <- ggarrange(plot_allocation_acq, plot_allocation_transf, ncol = 1, nrow = 2)

plot_gdp <- ggarrange(plot_gdp_acq, plot_gdp_transf, ncol = 1, nrow = 2)
plot_pop <- ggarrange(plot_pop_acq,plot_pop_transf, ncol = 1, nrow = 2)
plot_area <- ggarrange(plot_area_acq,plot_area_transf, ncol = 1, nrow = 2)


pdf("plot_core.pdf")
plot_core
dev.off()

pdf("plot_emission.pdf")
plot_emission
dev.off()

pdf("plot_allocation.pdf")
plot_allocation
dev.off()

pdf("plot_gdp.pdf")
plot_gdp
dev.off()

pdf("plot_pop.pdf")
plot_pop
dev.off()

pdf("plot_area.pdf")
plot_area
dev.off()


#############

estimates <- function(x,t){
  tmp <- x
  db_fitted <- as.data.frame(x$fitted.values)
  db_fitted$names <- rownames(db_fitted)
  db_residuals <- as.data.frame(x$residuals)
  db_residuals$names <- rownames(db_residuals)
  db_data <- as.data.frame(x$data)
  db_data$names <- rownames(db_data)
  db_all <- merge(db_fitted,db_residuals, by=c("names"))
  db_all <- merge(db_all,db_data[,c(1,2,ncol(db_data))], by=c("names"))
  db_all <- merge(db_all, db_model[,c("iso_o","iso_d",paste0("amount.",t))], by=c("iso_o","iso_d"))
  names(db_all)[6] <- "amount"
  db_all$year <- t
  return <- db_all
}

db_stats_2005 <- estimates(fit_2005, 2005)
db_stats_2006 <- estimates(fit_2006, 2006)
db_stats_2007 <- estimates(fit_2007, 2007)
db_stats_2008 <- estimates(fit_2008, 2008)
db_stats_2009 <- estimates(fit_2009, 2009)
db_stats_2010 <- estimates(fit_2010, 2010)
db_stats_2011 <- estimates(fit_2011, 2011)
db_stats_2012 <- estimates(fit_2012, 2012)
db_stats_2013 <- estimates(fit_2013, 2013)
db_stats_2014 <- estimates(fit_2014, 2014)
db_stats_2015 <- estimates(fit_2015, 2015)
db_stats_2016 <- estimates(fit_2016, 2016)
db_stats_2017 <- estimates(fit_2017, 2017)
db_stats_2018 <- estimates(fit_2018, 2018)
db_stats_2019 <- estimates(fit_2019, 2019)
db_stats_2020 <- estimates(fit_2020, 2020)

db_stats <- rbind(db_stats_2005,db_stats_2006,db_stats_2007,db_stats_2008,db_stats_2009,db_stats_2010,
                  db_stats_2011,db_stats_2012,db_stats_2013,db_stats_2014,db_stats_2015,db_stats_2016,
                  db_stats_2017,db_stats_2018,db_stats_2019,db_stats_2020)


write.csv(db_stats, "db_stats.csv")
write.csv(stats_db, "stats_db.csv")


db_deviance <- data.frame(fit_2005$deviance,fit_2006$deviance,fit_2007$deviance,fit_2008$deviance,
                          fit_2009$deviance,fit_2010$deviance,fit_2011$deviance,fit_2012$deviance,
                          fit_2013$deviance,fit_2014$deviance,fit_2015$deviance,fit_2016$deviance,
                          fit_2017$deviance,fit_2018$deviance,fit_2019$deviance,fit_2020$deviance)

write.csv(db_deviance, "db_deviance.csv")

xtable(stats_db[,c(1:9)])
xtable(stats_db[,c(1,10:17)])


