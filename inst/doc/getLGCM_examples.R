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
load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))

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
#  # Standardize time-invariant covariates (TICs)
#  ## ex1 and ex2 are standardized growth TICs in models
#  RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#  RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#  xstarts <- mean(baseT)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  Math_LGCM_BLS_f <- getLGCM(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#    intrinsic = TRUE, records = 1:9, growth_TIC = NULL, res_scale = 0.1
#  )
#  Math_LGCM_BLS_r <- getLGCM(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#    intrinsic = FALSE, records = 1:9, growth_TIC = NULL, res_scale = 0.1
#  )

## -----------------------------------------------------------------------------
getLRT(
  full = Math_LGCM_BLS_f@mxOutput, reduced = Math_LGCM_BLS_r@mxOutput, boot = FALSE, replications = NA
  )

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  paraBLS.TIC_LGCM.f <- c(
#    "alpha0", "alpha1", "alpha2", "alphag",
#    paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")),
#    "residuals",
#    paste0("beta1", c(0:2, "g")), paste0("beta2", c(0:2, "g")),
#    paste0("mux", 1:2), paste0("phi", c("11", "12", "22")),
#    "mueta0", "mueta1", "mueta2", "mu_knot"
#    )
#  Math_LGCM_TIC_BLS_f <- getLGCM(
#    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#    intrinsic = TRUE, records = 1:9, growth_TIC = c("ex1", "ex2"), res_scale = 0.1,
#    paramOut = TRUE, names = paraBLS.TIC_LGCM.f
#    )

## -----------------------------------------------------------------------------
Math_LGCM_TIC_BLS_f@Estimates
Figure1 <- getFigure(
  model = Math_LGCM_TIC_BLS_f@mxOutput, nClass = NULL, cluster_TIC = NULL, sub_Model = "LGCM",
  y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
  m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
  outcome = "Mathematics"
)
show(Figure1)

