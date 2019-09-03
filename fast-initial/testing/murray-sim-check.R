library(secr)
library(secrdesign)

nk <- function(x, k = 0.5, target= 20, grid, mask, noccasions = 5) {
  # x is lambda0; k is sigma in units of l_D
  # D = 1 because in lD units
  nrm <- Enrm(D = 10000, traps = grid, mask = mask,
              detectpar = list(lambda0 = x, sigma = k),
              noccasions = noccasions)
  nrm['Er'] - target
}

spacing <- 1
traps1 <- make.grid(nx = 8, ny = 8, detector = 'count', spacing = spacing)
mask <- make.mask(traps1, type = "trapbuffer", buffer = 4 * spacing)
noccasions <- 5
k <- 1 * spacing

# desired number of recaptures
r <- 50

# calculate lambda0 that gives you r recaptures in 1 occasion
ur <- try(uniroot(nk, interval = c(0.001, 100), target = r, k = k,
                  grid = traps1, mask = mask, noccasions = noccasions))
if(inherits(ur, 'try-error')){
  lam0 <- NA
}else{lam0 <- ur$root * noccasions}

set.seed(123)
# simulate capture history using that lambda0
CH <- sim.capthist(traps1, popn = list(D = 10000, buffer = 4 * spacing),
                   noccasions = 1, detectfn = "HHN",
                   detectpar = list(lambda0 = lam0, sigma = k))

nrow(CH) # n
sum(CH)-nrow(CH) # r

