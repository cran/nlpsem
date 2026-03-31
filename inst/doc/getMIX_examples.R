## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getMIX_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getMIX_examples.RData", package = "nlpsem"))

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
# ## gx1 and gx2 are standardized cluster TICs in models
# RMS_dat0$gx1 <- scale(RMS_dat0$INCOME)
# RMS_dat0$gx2 <- scale(RMS_dat0$EDU)
# xstarts <- mean(baseT)

## ----message = FALSE, eval = FALSE--------------------------------------------
# Math_BLS_LGCM1 <- getLGCM(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#   records = 1:9, tries = 10
#   )
# Math_BLS_LGCM2 <- getMIX(
#   dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM", y_var = "M",
#   t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE, tries = 10
# )
# set.seed(20191029)
# Math_BLS_LGCM3 <- getMIX(
#   dat = RMS_dat0, prop_starts = c(0.30, 0.40, 0.30), sub_Model = "LGCM", y_var = "M",
#   t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE, tries = 10
# )

## -----------------------------------------------------------------------------
Figure1 <- getFigure(
  model = Math_BLS_LGCM1@mxOutput, nClass = NULL, cluster_TIC = NULL, sub_Model = "LGCM",
  y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
  m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
  outcome = "Mathematics"
)
show(Figure1)
Figure2 <- getFigure(
  model = Math_BLS_LGCM2@mxOutput, nClass = 2, cluster_TIC = NULL, sub_Model = "LGCM",
  y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
  m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
  outcome = "Mathematics"
)
show(Figure2)
Figure3 <- getFigure(
  model = Math_BLS_LGCM3@mxOutput, nClass = 3, cluster_TIC = NULL, sub_Model = "LGCM",
  y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
  m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
  outcome = "Mathematics"
)
show(Figure3)
getSummary(model_list = list(Math_BLS_LGCM1@mxOutput, Math_BLS_LGCM2@mxOutput, Math_BLS_LGCM3@mxOutput),
           HetModels = TRUE)

## ----message = FALSE, eval = FALSE--------------------------------------------
# set.seed(20191029)
# RM_BLS_PLGCM3 <- getMIX(
#   dat = RMS_dat0, prop_starts = c(0.30, 0.40, 0.30), sub_Model = "MGM",
#   cluster_TIC = c("gx1", "gx2"), t_var = c("T", "T"), y_var = c("R", "M"),
#   curveFun = "BLS", intrinsic = FALSE, records = list(1:9, 1:9),
#   y_model = "LGCM", tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Figure4 <- getFigure(
  model = RM_BLS_PLGCM3@mxOutput, nClass = 3, cluster_TIC = c("gx1", "gx2"), 
  sub_Model = "MGM", y_var = c("R", "M"), curveFun = "BLS", y_model = "LGCM",
  t_var = c("T", "T"), records = list(1:9, 1:9), m_var = NULL, x_var = NULL, 
  x_type = NULL, xstarts = xstarts, xlab = "Month", 
  outcome = c("Reading", "Mathematics")
)
show(Figure4)

