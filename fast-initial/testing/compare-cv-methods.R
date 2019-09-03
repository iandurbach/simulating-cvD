library(dplyr)
library(purrr)
library(secr)
library(secrdesign)

spacing <- 1000
traps1 <- make.grid(nx = 4, ny = 4, detector = 'count', spacing = spacing)
mask <- make.mask(traps1, type = "trapbuffer", buffer = 4 * spacing)
noccasions <- 1
sigma <- 1 * spacing
lam0 <- 0.08717
D <- 0.06

### manually run once

# simulate capture history using that lambda0
CH <- sim.capthist(traps1, popn = list(D = D, buffer = 4 * spacing),
                   noccasions = noccasions, detectfn = "HHN",
                   detectpar = list(lambda0 = lam0, sigma = sigma))

nrow(CH) # n
sum(CH)-nrow(CH) # r

# fit secr model
fit <- secr.fit(CH, mask = mask, detectfn = "HHN", trace = F)
out <- predict(fit) %>% mutate(par = c("D", "lambda0", "sigma"),
                               cv = SE.estimate / estimate,
                               n = nrow(CH),
                               r = sum(CH)-nrow(CH),
                               rotCV = 1/sqrt(pmin(n,r)))
out

### manually run many times

# first bundle stuff above into function

sim_and_fit <- function(traps, D, spacing, noccasions, lam0, sigma, mask){

  # simulate capture history using that lambda0
  CH <- sim.capthist(traps, popn = list(D = D, buffer = 4 * spacing),
                     noccasions = noccasions, detectfn = "HHN",
                     detectpar = list(lambda0 = lam0, sigma = sigma))

  nrow(CH) # n
  sum(CH)-nrow(CH) # r

  # fit secr model
  fit <- secr.fit(CH, mask = mask, detectfn = "HHN", trace = F)
  out <- predict(fit) %>% mutate(par = c("D", "lambda0", "sigma"),
                                 cv = SE.estimate / estimate,
                                 n = nrow(CH),
                                 r = sum(CH)-nrow(CH),
                                 rotCV = 1/sqrt(pmin(n,r))) %>%
    filter(par == "D")

  return(out)
}

# run
m0 <- 10 %>%
  rerun(sim_and_fit(traps = traps1, D = D, spacing = spacing, noccasions = noccasions,
                    lam0 = lam0, sigma = sigma, mask = mask)) %>%
  bind_rows()

res0 <- m0 %>% group_by(par) %>%
  summarize(m0_n_mean = mean(n), m0_r_mean = mean(r),
            m0_est_mean = mean(estimate), m0_se.est_mean = mean(SE.estimate),
            m0_cv_mean = mean(cv), m0_cv_median = median(cv),
            m0_rotCV_mean = mean(rotCV)) %>% data.frame()
res0

### CV(D) the rule-of-thumb way -- 1/sqrt(min(n,r))

mydetectors <- list(traps1)

scen <- make.scenarios(trapsindex = 1, D = D,
                       sigma = sigma, nrepeats = 1, detectfn = "HHN",
                       noccasions = noccasions, lambda0 = lam0)


m1 <- scenarioSummary(scen, mydetectors)

res1 <- m1 %>% mutate(par = "D") %>% select(par, m1_En = En, m1_Er = Er, m1_rotCV = rotRSE)

res1

### CV(D) the intermediate way -- gradient of likelihood at assumed par values

m2 <- run.scenarios(nrepl = 50, trapset = mydetectors, scenarios = scen,
                      seed = 345, fit = TRUE,
                      fit.args = list(method = "none"),
                      extract = derived)

summary(m2, fields = c("n", "mean", "median", "se"))

res2 <- data.frame(par = "D",
                   m2_cv_mean = summary(m2, fields = c("n", "mean", "median", "se"))$OUTPUT$`1`[7,2],
                   m2_cv_median = summary(m2, fields = c("n", "mean", "median", "se"))$OUTPUT$`1`[7,3],
                   m2_cv_se = summary(m2, fields = c("n", "mean", "median", "se"))$OUTPUT$`1`[7,4])

res2

### CV(D) the hard way -- full simulation approach

m3 <- run.scenarios(nrepl = 20, trapset = mydetectors, scenarios = scen,
                      seed = 345, fit = TRUE,
                      fit.function = c("secr.fit"),
                      fit.args = list(detectfn = "HHN"),
                      extract = derived)

summary(m3, fields = c("n", "mean", "median", "se"))

res3 <- data.frame(par = "D",
                   m3_cv_mean = summary(m3, fields = c("n", "mean", "median", "se"))$OUTPUT$`1`[7,2],
                   m3_cv_median = summary(m3, fields = c("n", "mean", "median", "se"))$OUTPUT$`1`[7,3],
                   m3_cv_se = summary(m3, fields = c("n", "mean", "median", "se"))$OUTPUT$`1`[7,4])
res3

### put all methods together
pars <- data.frame(D = D, spacing = spacing, noccasions = noccasions, lam0 = lam0, sigma = sigma)
res_all <- pars %>% cbind(res0) %>% left_join(res1, by = "par") %>% left_join(res2, by = "par") %>% left_join(res3, by = "par")
res_all
