## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(OpenMx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ---- message = FALSE---------------------------------------------------------
load(system.file("extdata", "getTVCmodel_examples.RData", package = "nlpsem"))

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  # Load ECLS-K (2011) data
#  data("RMS_dat")
#  RMS_dat0 <- RMS_dat
#  # Re-baseline the data so that the estimated initial status is for the
#  # starting point of the study
#  baseT <- RMS_dat0$T1
#  RMS_dat0$T1 <- (RMS_dat0$T1 - baseT)/12
#  RMS_dat0$T2 <- (RMS_dat0$T2 - baseT)/12
#  RMS_dat0$T3 <- (RMS_dat0$T3 - baseT)/12
#  RMS_dat0$T4 <- (RMS_dat0$T4 - baseT)/12
#  RMS_dat0$T5 <- (RMS_dat0$T5 - baseT)/12
#  RMS_dat0$T6 <- (RMS_dat0$T6 - baseT)/12
#  RMS_dat0$T7 <- (RMS_dat0$T7 - baseT)/12
#  RMS_dat0$T8 <- (RMS_dat0$T8 - baseT)/12
#  RMS_dat0$T9 <- (RMS_dat0$T9 - baseT)/12
#  # Standardize time-invariant covariates (TICs)
#  ## ex1 is standardized growth TIC in models
#  RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#  # Standardize time-varying covariate (TVC)
#  BL_mean <- mean(RMS_dat0[, "R1"])
#  BL_var <- var(RMS_dat0[, "R1"])
#  RMS_dat0$Rs1 <- (RMS_dat0$R1 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs2 <- (RMS_dat0$R2 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs3 <- (RMS_dat0$R3 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs4 <- (RMS_dat0$R4 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs5 <- (RMS_dat0$R5 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs6 <- (RMS_dat0$R6 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs7 <- (RMS_dat0$R7 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs8 <- (RMS_dat0$R8 - BL_mean)/sqrt(BL_var)
#  RMS_dat0$Rs9 <- (RMS_dat0$R9 - BL_mean)/sqrt(BL_var)
#  xstarts <- mean(baseT)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  set.seed(20191029)
#  Math_TVC_BLS_f <- getTVCmodel(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE,
#    records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 0, growth_TIC = "ex1",
#    res_scale = 0.1, tries = 10
#    )
#  paraBLS_TVC.f <- c(
#    "Y_alpha0", "Y_alpha1", "Y_alpha2", "Y_alphag",
#    paste0("Y_psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")),
#    "Y_residuals", "X_mueta0", "X_mueta1", paste0("X_psi", c("00", "01", "11")),
#    paste0("X_rel_rate", 2:8), paste0("X_abs_rate", 1:8), "X_residuals",
#    paste0("betaTIC", c(0:2, "g")), paste0("betaTVC", c(0:2, "g")), "muTIC", "phiTIC",
#    "Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_mu_knot", "covBL", "kappa", "Cov_XYres"
#    )
#  set.seed(20191029)
#  Math_TVCslp_BLS_f <- getTVCmodel(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE,
#    records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 1, growth_TIC = "ex1",
#    res_scale = c(0.1, 0.1), res_cor = 0.3, tries = 10, paramOut = TRUE,
#    names = paraBLS_TVC.f
#    )

## -----------------------------------------------------------------------------
getEstimateStats(est_in = Math_TVCslp_BLS_f[[2]], CI_type = "Wald")
Figure1 <- getFigure(
  model = Math_TVC_BLS_f, sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
print(Figure1)
Figure2 <- getFigure(
  model = Math_TVCslp_BLS_f[[1]], sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
print(Figure2)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  paraBLS_TVC.r <- c(
#    "Y_alpha0", "Y_alpha1", "Y_alpha2", "Y_knot",
#    paste0("Y_psi", c("00", "01", "02", "11", "12", "22")), "Y_residuals",
#    "X_mueta0", "X_mueta1", paste0("X_psi", c("00", "01", "11")),
#    paste0("X_rel_rate", 2:8), paste0("X_abs_rate", 1:8), "X_residuals",
#    paste0("betaTIC", 0:2), paste0("betaTVC", 0:2), "muTIC", "phiTIC",
#    "Y_mueta0", "Y_mueta1", "Y_mueta2", "covBL", "kappa", "Cov_XYres"
#    )
#  set.seed(20191029)
#  Math_TVCslp_BLS_r <- getTVCmodel(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#    records = 1:9, y_model = "LGCM", TVC = "R", decompose = 1, growth_TIC = "ex1",
#    res_scale = c(0.1, 0.1), res_cor = 0.3, tries = 10,  paramOut = TRUE,
#    names = paraBLS_TVC.r)
#  set.seed(20191029)
#  Math_TVCchg_BLS_r <- getTVCmodel(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#    records = 1:9, y_model = "LGCM", TVC = "R", decompose = 2, growth_TIC = "ex1",
#    res_scale = c(0.1, 0.1), res_cor = 0.3, tries = 10,  paramOut = TRUE,
#    names = paraBLS_TVC.r)
#  set.seed(20191029)
#  Math_TVCchgBL_BLS_r <- getTVCmodel(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#    records = 1:9, y_model = "LGCM", TVC = "R", decompose = 3, growth_TIC = "ex1",
#    res_scale = c(0.1, 0.1), res_cor = 0.3, tries = 10,  paramOut = TRUE,
#    names = paraBLS_TVC.r)

## -----------------------------------------------------------------------------
Figure3 <- getFigure(
  model = Math_TVCslp_BLS_r[[1]], sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
print(Figure3)
Figure4 <- getFigure(
  model = Math_TVCchg_BLS_r[[1]], sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
print(Figure4)
Figure5 <- getFigure(
  model = Math_TVCchgBL_BLS_r[[1]], sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
print(Figure5)

