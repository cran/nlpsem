## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getMGroup_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getMGroup_examples.RData", package = "nlpsem"))

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
# set.seed(20191029)
# MGroup_Math_BLS_LGCM_f <-  getMGroup(
#   dat = RMS_dat0, grp_var = "SEX", sub_Model = "LGCM", y_var = "M", t_var = "T",
#   records = 1:9, curveFun = "BLS", intrinsic = TRUE, tries = 20
# )

## -----------------------------------------------------------------------------
Figure1 <- getFigure(
  model = MGroup_Math_BLS_LGCM_f@mxOutput, nClass = 2, cluster_TIC = NULL, grp_var = "SEX",
  sub_Model = "LGCM", y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T",
  records = 1:9, m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts,
  xlab = "Month", outcome = "Mathematics"
)
show(Figure1)

## ----message = FALSE, eval = FALSE--------------------------------------------
# set.seed(20191029)
# MGroup_Read_EXP_LGCM_TIC_r <- getMGroup(
#   dat = RMS_dat0, grp_var = "SEX", sub_Model = "LGCM", y_var = "R", t_var = "T",
#   records = 1:9, curveFun = "EXP", intrinsic = FALSE,
#   growth_TIC = c("ex1", "ex2"), paramOut = TRUE
# )

## ----eval = FALSE-------------------------------------------------------------
# MGroup_Read_EXP_LGCM_TIC_r@Estimates
# Figure2 <- getFigure(
#   model = MGroup_Read_EXP_LGCM_TIC_r@mxOutput, nClass = 2, cluster_TIC = NULL, grp_var = "SEX",
#   sub_Model = "LGCM", y_var = "R", curveFun = "EXP", y_model = "LGCM", t_var = "T",
#   records = 1:9, m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts,
#   xlab = "Month", outcome = "Reading"
# )
# show(Figure2)

