## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ---- message = FALSE---------------------------------------------------------
load(system.file("extdata", "getMediation_examples.RData", package = "nlpsem"))

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  # Load ECLS-K (2011) data
#  data("RMS_dat")
#  RMS_dat0 <- RMS_dat
#  # Re-baseline the data so that the estimated initial status is for the
#  # starting point of the study
#  baseT <- RMS_dat0$T1
#  RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#  RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#  RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#  RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#  RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#  RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#  RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#  RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#  RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#  RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#  xstarts <- mean(baseT)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  paraMed2_BLS <- c(
#    "muX", "phi11", "alphaM1", "alphaMr", "alphaM2", "mugM",
#    paste0("psi", c("M1M1", "M1Mr", "M1M2", "MrMr", "MrM2", "M2M2"), "_r"),
#    "alphaY1", "alphaYr", "alphaY2", "mugY",
#    paste0("psi", c("Y1Y1", "Y1Yr", "Y1Y2", "YrYr", "YrY2", "Y2Y2"), "_r"),
#    paste0("beta", rep(c("M", "Y"), each = 3), rep(c(1, "r", 2), 2)),
#    paste0("beta", c("M1Y1", "M1Yr", "M1Y2", "MrYr", "MrY2", "M2Y2")),
#    "muetaM1", "muetaMr", "muetaM2", "muetaY1", "muetaYr", "muetaY2",
#    paste0("Mediator", c("11", "1r", "12", "rr", "r2", "22")),
#    paste0("total", c("1", "r", "2")),
#    "residualsM", "residualsY", "residualsYM"
#    )
#  Med2_LGCM_BLS <- getMediation(
#    dat = RMS_dat0, t_var = rep("T", 2), y_var = "M", m_var = "R",
#    x_type = "baseline", x_var = "ex1", curveFun = "bilinear spline",
#    records = list(1:9, 1:9), res_scale = c(0.1, 0.1), res_cor = 0.3,
#    paramOut = TRUE, names = paraMed2_BLS
#    )

## -----------------------------------------------------------------------------
Med2_LGCM_BLS@Estimates

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  paraMed3_BLS <- c(
#    "muetaX1", "muetaXr", "muetaX2", "mugX",
#    paste0("psi", c("X1X1", "X1Xr", "X1X2", "XrXr", "XrX2", "X2X2")),
#    "alphaM1", "alphaMr", "alphaM2", "mugM",
#    paste0("psi", c("M1M1", "M1Mr", "M1M2", "MrMr", "MrM2", "M2M2"), "_r"),
#    "alphaY1", "alphaYr", "alphaY2", "mugY",
#    paste0("psi", c("Y1Y1", "Y1Yr", "Y1Y2", "YrYr", "YrY2", "Y2Y2"), "_r"),
#    paste0("beta", c("X1Y1", "X1Yr", "X1Y2", "XrYr", "XrY2", "X2Y2",
#                     "X1M1", "X1Mr", "X1M2", "XrMr", "XrM2", "X2M2",
#                     "M1Y1", "M1Yr", "M1Y2", "MrYr", "MrY2", "M2Y2")),
#    "muetaM1", "muetaMr", "muetaM2", "muetaY1", "muetaYr", "muetaY2",
#    paste0("mediator", c("111", "11r", "112", "1rr", "1r2", "122", "rr2", "r22", "rrr", "222")),
#    paste0("total", c("11", "1r", "12", "rr", "r2", "22")),
#    "residualsX", "residualsM", "residualsY", "residualsMX", "residualsYX", "residualsYM"
#    )
#  set.seed(20191029)
#  Med3_LGCM_BLS <- getMediation(
#    dat = RMS_dat0, t_var = rep("T", 3), y_var = "S", m_var = "M", x_type = "longitudinal",
#    x_var = "R", curveFun = "bilinear spline", records = list(2:9, 1:9, 1:9),
#    res_scale = c(0.1, 0.1, 0.1),  res_cor = c(0.3, 0.3), tries = 10, paramOut = TRUE,
#    names = paraMed3_BLS
#    )

## -----------------------------------------------------------------------------
Med3_LGCM_BLS@Estimates

