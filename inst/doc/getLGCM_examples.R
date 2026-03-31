## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))

## ----message = FALSE, eval = FALSE--------------------------------------------
# # Load ECLS-K (2011) data
# data("RMS_dat")
# RMS_dat0 <- RMS_dat
# # Re-baseline the data so that the estimated initial status is for the
# # starting point of the study
# baseT <- RMS_dat0$T1
# RMS_dat0$T1 <- RMS_dat0$T1 - baseT
# RMS_dat0$T2 <- RMS_dat0$T2 - baseT
# RMS_dat0$T3 <- RMS_dat0$T3 - baseT
# RMS_dat0$T4 <- RMS_dat0$T4 - baseT
# RMS_dat0$T5 <- RMS_dat0$T5 - baseT
# RMS_dat0$T6 <- RMS_dat0$T6 - baseT
# RMS_dat0$T7 <- RMS_dat0$T7 - baseT
# RMS_dat0$T8 <- RMS_dat0$T8 - baseT
# RMS_dat0$T9 <- RMS_dat0$T9 - baseT
# # Standardize time-invariant covariates (TICs)
# ## ex1 and ex2 are standardized growth TICs in models
# RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
# RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
# xstarts <- mean(baseT)

## ----message = FALSE, eval = FALSE--------------------------------------------
# Math_LGCM_BLS_f <- getLGCM(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#   intrinsic = TRUE, records = 1:9, growth_TIC = NULL, tries = 10
# )
# Math_LGCM_BLS_r <- getLGCM(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#   intrinsic = FALSE, records = 1:9, growth_TIC = NULL, tries = 10
# )

## -----------------------------------------------------------------------------
getLRT(
  full = Math_LGCM_BLS_f@mxOutput, reduced = Math_LGCM_BLS_r@mxOutput, boot = FALSE, rep = NA
  )

## ----message = FALSE, eval = FALSE--------------------------------------------
# Math_LGCM_TIC_BLS_f <- getLGCM(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#   intrinsic = TRUE, records = 1:9, growth_TIC = c("ex1", "ex2"),
#   tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Math_LGCM_TIC_BLS_f@Estimates
Figure1 <- getFigure(
  model = Math_LGCM_TIC_BLS_f@mxOutput, nClass = NULL, cluster_TIC = NULL, sub_Model = "LGCM",
  y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
  m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
  outcome = "Mathematics"
)
show(Figure1)

