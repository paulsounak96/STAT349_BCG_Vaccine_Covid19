library(geosphere)
library(regress)
library(dplyr)

ttestdata = cumweeklydeaths %>% 
  group_by(Country, ISO3, population_million, HDI_2018, log_pop_density,
           ages_65_up, BCG_policy, BCG_mean_coverage, range_age_BCG, BCG_index,
           urban_fraction, lat, long) %>%
  summarise(deaths = Deaths[32], response = log_deaths_million[32])
ttestdata = as.data.frame(ttestdata[!is.na(ttestdata$response), ])

spatial_cov = matrix(1, 142, 142)
for (i in 1:142) {
  for (j in 1:142) {
    spatial_cov[i,j]=exp(-distHaversine(c(ttestdata$long[i],ttestdata$lat[i]),
                                        c(ttestdata$long[j],ttestdata$lat[j]))/20037508)
  }
}

x0 = model.matrix(~ ttestdata$HDI_2018 + ttestdata$log_pop_density +
                    ttestdata$ages_65_up + ttestdata$urban_fraction)

loglikelihood = c()

for (rho in seq(0.0001, 1, 0.001)) {
  m = spatial_cov^rho
  mod_notemp0 = regress(response ~ HDI_2018 + log_pop_density + ages_65_up +
                          urban_fraction, ~m, start = mod_notemp0$sigma,
                        data = ttestdata)
  loglikelihood = c(loglikelihood, mod_notemp0$llik)
}
plot(seq(0.0001, 1, 0.001), loglikelihood)

#-----------



###########################rho goes to 0##############################

spatial_cov1 = matrix(1, 142, 142)
for (i in 1:142) {
  for (j in 1:142) {
    spatial_cov1[i,j]=-distHaversine(c(ttestdata$long[i],ttestdata$lat[i]),
                                        c(ttestdata$long[j],ttestdata$lat[j]))/20037508
  }
}

x0 = model.matrix(~ ttestdata$HDI_2018 + ttestdata$log_pop_density +
                    ttestdata$ages_65_up + ttestdata$urban_fraction)

loglikelihood1 = c()

for (rho in seq(0.0005, 0.1, 0.0005)) {
  m = rho*spatial_cov1
  mod_notemp0 = regress(response ~ HDI_2018 + log_pop_density + ages_65_up +
                          urban_fraction, ~m, 
                        start = mod_notemp0$sigma, data = ttestdata)
  loglikelihood1 = c(loglikelihood1, mod_notemp0$llik)
}
plot(seq(0.0005, 0.1, 0.0005), loglikelihood1)

###########################Spatial Only##############################

spatial_cov1 = matrix(1, 142, 142)
for (i in 1:142) {
  for (j in 1:142) {
    spatial_cov1[i,j]=-distHaversine(c(ttestdata$long[i],ttestdata$lat[i]),
                                     c(ttestdata$long[j],ttestdata$lat[j]))/20037508
  }
}

x0 = model.matrix(~ ttestdata$HDI_2018 + ttestdata$log_pop_density +
                    ttestdata$ages_65_up + ttestdata$urban_fraction)

mod_ttest0 = regress(response ~ HDI_2018 + log_pop_density + ages_65_up +
                          urban_fraction, ~spatial_cov1, 
                        start = mod_notemp0$sigma, data = ttestdata)

mod_ttest1a = regress(response ~ HDI_2018 + log_pop_density + ages_65_up +
                        urban_fraction + BCG_index, ~spatial_cov1, kernel = x0,
                      start = mod_notemp0$sigma, data = ttestdata)

2*(mod_ttest1a$llik - mod_ttest0$llik)
pchisq(2*(mod_ttest1a$llik - mod_ttest0$llik), 1, lower.tail = FALSE)

mod_ttest1b = regress(response ~ HDI_2018 + log_pop_density + ages_65_up +
                         BCG_index + urban_fraction, ~spatial_cov1,
                       start = mod_ttest0$sigma, data = ttestdata)
summary(mod_ttest1b)

#######################Spatiotemporal Models############################

ctry = c(as.character(sample(countryinfo$Country[countryinfo$BCG_policy == "current"], 70)),
         as.character(countryinfo$Country[countryinfo$BCG_policy != "current"]))

countryinfo = countryinfo[is.element(countryinfo$Country, ctry), ]
cumweeklydeaths = cumweeklydeaths[is.element(cumweeklydeaths$Country, ctry), ]

cumweeklydeaths = inner_join(data.frame(id = 1:99, Country= countryinfo$Country),
                           cumweeklydeaths, by = "Country")
cumweeklydeaths = cumweeklydeaths[cumweeklydeaths$Week <= 15, ]

sm = matrix(1, 99, 99)
exp_sm = matrix(1, 99, 99)
spatial_cov = matrix(1, 1485, 1485)
exp_spatial_cov = matrix(1, 1485, 1485)

for (i in 1:99) {
  for (j in 1:99) {
    sm[i,j] = -distHaversine(c(countryinfo$long[i],countryinfo$lat[i]),
                                 c(countryinfo$long[j],countryinfo$lat[j]))/20037508
  }
}
for (i in 1:99) {
  for (j in 1:99) {
    exp_sm[i,j] = exp(-distHaversine(c(countryinfo$long[i],countryinfo$lat[i]),
                             c(countryinfo$long[j],countryinfo$lat[j]))/20037508)
  }
}
for (i in 1:1485) {
  for (j in 1:1485) {
    spatial_cov[i,j] = sm[cumweeklydeaths$id[i], cumweeklydeaths$id[j]]
  }
}
for (i in 1:1485) {
  for (j in 1:1485) {
    exp_spatial_cov[i,j] = exp_sm[cumweeklydeaths$id[i], cumweeklydeaths$id[j]]
  }
}

temporal_cov = -abs(outer(cumweeklydeaths$Week, cumweeklydeaths$Week, "-"))
st_cov = exp_spatial_cov*exp(temporal_cov)

x0 = model.matrix(~ cumweeklydeaths$Week + cumweeklydeaths$HDI_2018 +
                    cumweeklydeaths$log_pop_density + 
                    cumweeklydeaths$ages_65_up + cumweeklydeaths$urban_fraction)

############################Don't use model###################################

mod_notemp0 = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                        ages_65_up + urban_fraction, ~ st_cov, tol = 0.001,
                      pos = c(1,1), data = cumweeklydeaths)

mod_notemp1a = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ st_cov,
                       kernel = x0, pos = c(1,1), start = mod_notemp0$sigma,
                       tol = 0.001, data = cumweeklydeaths)


#########################More general model###################################

mod_notemp2 = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                        ages_65_up + urban_fraction, ~ spatial_cov + temporal_cov, tol = 0.001,
                      data = cumweeklydeaths)

mod_notemp3a = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ spatial_cov + temporal_cov,
                       kernel = x0, start = mod_notemp2$sigma,
                       tol = 0.001, data = cumweeklydeaths)

2*(mod_notemp3a$llik - mod_notemp2$llik)
pchisq(2*(mod_notemp3a$llik - mod_notemp2$llik), 1, lower.tail = FALSE)

mod_notemp3b = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ spatial_cov + temporal_cov,
                       start = mod_notemp2$sigma,
                       tol = 0.001, data = cumweeklydeaths)
summary(mod_notemp3b)

#########################Most general model###################################

mod_notemp4 = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                        ages_65_up + urban_fraction, ~ spatial_cov + temporal_cov + st_cov,
                      tol = 0.001, data = cumweeklydeaths)

mod_notemp5a = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ spatial_cov + temporal_cov + st_cov,
                       kernel = x0, start = mod_notemp4$sigma,
                       tol = 0.001, data = cumweeklydeaths)

2*(mod_notemp5a$llik - mod_notemp4$llik)
pchisq(2*(mod_notemp5a$llik - mod_notemp4$llik), 1, lower.tail = FALSE)

mod_notemp5b = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ spatial_cov + temporal_cov + st_cov,
                       start = mod_notemp4$sigma,
                       tol = 0.001, data = cumweeklydeaths)
summary(mod_notemp5b)

#########################Final model###################################

mod_notemp6 = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                        ages_65_up + urban_fraction, ~ spatial_cov + st_cov,
                      tol = 0.001, data = cumweeklydeaths)

mod_notemp7a = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ spatial_cov + st_cov,
                       kernel = x0, start = mod_notemp6$sigma, tol = 0.001, data = cumweeklydeaths)

2*(mod_notemp7a$llik - mod_notemp6$llik)
pchisq(2*(mod_notemp7a$llik - mod_notemp6$llik), 1, lower.tail = FALSE)

mod_notemp7b = regress(log_deaths_million ~ Week + HDI_2018 + log_pop_density +
                         BCG_index + ages_65_up + urban_fraction, ~ spatial_cov + st_cov,
                       start = mod_notemp6$sigma, tol = 0.001, data = cumweeklydeaths)
summary(mod_notemp7b)