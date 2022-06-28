"
Author/maintainer: Can Zhou [eidotog@gmail.com]
Date: Oct 4 2021
Version: 0.2
"

dat <- read.csv("../data/_tagBFT.csv")

dat.c <- dat[!is.na(as.numeric(dat$ReLenCM)) & !is.na(as.numeric(dat$RcLenCM)) & !is.na(as.numeric(dat$RelYear)) & !is.na(as.numeric(dat$RecYear)), ]

west <- which(dat.c$ReLonX < -45)
east <- which(dat.c$ReLonX > -45)

dat.c$reg <- 0
dat.c$reg[east] <- 1

t1 <- as.Date(paste(dat.c$RelYear, dat.c$RelMon, dat.c$RelDay, sep = "-"), "%Y-%m-%d")
ym <- zoo::as.yearmon(t1)
t1_next_month <- zoo::as.Date(ym, frac = 1) + 1
t2 <- as.Date(paste(dat.c$RecYear, dat.c$RecMon, dat.c$RecDay, sep = "-"), "%Y-%m-%d")
ym2 <- zoo::as.yearmon(t2)
t2_next_month <- zoo::as.Date(ym2, frac = 1) + 1

d1 <- as.numeric(difftime(t1_next_month, t1, units = "days"))
d2 <- as.numeric(difftime(t2_next_month, t2, units = "days"))

dat.c$d1 <- 1 - d1 / 30
dat.c$d2 <- 1 - d2 / 30

day1 <- as.Date("1950-1-1", "%Y-%m-%d")
d1.d <- as.integer(difftime(t1, day1, units = "days"))
d2.d <- as.integer(difftime(t2, day1, units = "days"))

dat.c$d1.d <- d1.d
dat.c$d2.d <- d2.d

year.min <- min(as.numeric(dat.c$RelYear)) - 1 # 1960

nao.d <- read.csv("../data/norm.daily.nao.index.b500101.current.ascii", sep = "", header = FALSE)
nao.d$Date <- 1:(dim(nao.d)[1]) # 1950-01-01 => Day 1

ao <- read.csv("../data/monthly.ao.index.b50.current.ascii", sep = "", header = FALSE)
ao$Date <- (ao$V1 - year.min) * 12 + ao$V2

pna <- read.csv("../data/norm.pna.monthly.b5001.current.ascii", sep = "", header = FALSE)
pna$Date <- (pna$V1 - year.min) * 12 + pna$V2

nao <- read.csv("../data/norm.nao.monthly.b5001.current.ascii", sep = "", header = FALSE)
nao$Date <- (nao$V1 - year.min) * 12 + nao$V2

amo_ <- read.csv("../data/amon", sep = "", header = FALSE)
tmp1 <- numeric()
for (i in 1:dim(amo_)[1]) {
  tmp1 <- c(tmp1, as.numeric(amo_[i, 2:13]))
}
tmp0 <- numeric()
for (i in 1:dim(amo_)[1]) {
  tmp0 <- c(tmp0, rep(amo_$V1[i], 12))
}
amo <- data.frame(V1 = tmp0, V2 = 1:12, V3 = tmp1)
amo$Date <- (amo$V1 - year.min) * 12 + amo$V2
amo <- amo[!is.na(amo$V3), ]

fit0_ <- lm(ao$V3 ~ nao$V3)
ao$res <- resid(fit0_)
fit1_ <- lm(pna$V3 ~ nao$V3 + ao$V3)
pna$res <- resid(fit1_)

ao.c <- ao[which(ao$Date == 0):dim(ao)[1], ]
pna.c <- pna[which(pna$Date == 0):dim(pna)[1], ]
nao.c <- nao[which(nao$Date == 0):dim(nao)[1], ]
amo.c <- amo[which(amo$Date == 0):dim(amo)[1], ]

dat.c$REL_Date <- (as.numeric(dat.c$RelYear) - year.min) * 12 + as.numeric(dat.c$RelMon)
dat.c$REC_Date <- (as.numeric(dat.c$RecYear) - year.min) * 12 + as.numeric(dat.c$RecMon)

n_this_month <- as.Date(paste(nao$V1, nao$V2, "01", sep = "-"), "%Y-%m-%d")
n_ym <- zoo::as.yearmon(n_this_month)
n_next_month <- zoo::as.Date(n_ym, frac = 1) + 1
n_days <- as.numeric(difftime(n_next_month, n_this_month, units = "days"))
n_days.c <- n_days[which(nao$Date == 0):dim(nao)[1]]

ex <- which(dat.c$d2.d < dat.c$d1.d)
same_day <- which(dat.c$d1.d == dat.c$d2.d)
no_day <- which(is.na(dat.c$d1) | is.na(dat.c$d2))
too_long <- which(dat.c$d2.d - dat.c$d1.d > 365 * 20)
too_small <- which(dat.c$RcLenCM < 10)
exx <- union(union(union(union(ex, same_day), no_day), too_long), too_small)
datx_ <- dat.c[-exx, ]

fl <- which(datx_$RcLenType %in% c("FL", "LJF") & datx_$ReLenType %in% c("FL", "LJF"))
tl <- which(datx_$RcLenType %in% c("TLE") | datx_$ReLenType %in% c("TLE"))

tmpx <- datx_[fl, ]

## Climate index
nkts <- 21

gam_design <- mgcv::gam(V3 ~ s(V3, bs = "cs", k = nkts), data = nao, fit = FALSE)
gam_design$smooth[[1]]$xp
S <- gam_design$smooth[[1]]$S[[1]]
Sdim <- dim(S)[1]

nao.p <- c(0, seq(min(nao$V3, na.rm = TRUE), max(nao$V3, na.rm = TRUE), length.out = 100))
nao.p <- sort(nao.p)
P_design <- mgcv::PredictMat(gam_design$smooth[[1]], data = data.frame(V3 = nao.p))
X <- mgcv::PredictMat(gam_design$smooth[[1]], data = data.frame(V3 = nao.c$V3))
p_dm_nao <- mgcv::PredictMat(gam_design$smooth[[1]], data = data.frame(V3 = nao$V3))

if (MULTI) {
  ao$V3 <- ao$res
  ao.c$V3 <- ao.c$res
}

gam_design_a <- mgcv::gam(V3 ~ s(V3, bs = "cs", k = nkts), data = ao, fit = FALSE)
gam_design_a$smooth[[1]]$xp
S_a <- gam_design_a$smooth[[1]]$S[[1]]
Sdim_a <- dim(S_a)[1]

ao.p <- c(0, seq(min(ao$V3, na.rm = TRUE), max(ao$V3, na.rm = TRUE), length.out = 100))
ao.p <- sort(ao.p)
P_design_a <- mgcv::PredictMat(gam_design_a$smooth[[1]], data = data.frame(V3 = ao.p))
X_a <- mgcv::PredictMat(gam_design_a$smooth[[1]], data = data.frame(V3 = ao.c$V3))
p_dm_ao <- mgcv::PredictMat(gam_design_a$smooth[[1]], data = data.frame(V3 = ao$V3))

if (MULTI) {
  pna$V3 <- pna$res
  pna.c$V3 <- pna.c$res
}

gam_design_p <- mgcv::gam(V3 ~ s(V3, bs = "cs", k = nkts), data = pna, fit = FALSE)
gam_design_p$smooth[[1]]$xp
S_p <- gam_design_p$smooth[[1]]$S[[1]]
Sdim_p <- dim(S_p)[1]

pna.p <- c(0, seq(min(pna$V3, na.rm = TRUE), max(pna$V3, na.rm = TRUE), length.out = 100))
pna.p <- sort(pna.p)
P_design_p <- mgcv::PredictMat(gam_design_p$smooth[[1]], data = data.frame(V3 = pna.p))
X_p <- mgcv::PredictMat(gam_design_p$smooth[[1]], data = data.frame(V3 = pna.c$V3))
p_dm_pna <- mgcv::PredictMat(gam_design_p$smooth[[1]], data = data.frame(V3 = pna$V3))