## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getTVCmodel_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getTVCmodel_examples.RData", package = "nlpsem"))

## ----message = FALSE, eval = FALSE--------------------------------------------
# # Load ECLS-K (2011) data
# data("RMS_dat")
# RMS_dat0 <- RMS_dat
# # Re-baseline the data so that the estimated initial status is for the
# # starting point of the study
# baseT <- RMS_dat0$T1
# RMS_dat0$T1 <- (RMS_dat0$T1 - baseT)/12
# RMS_dat0$T2 <- (RMS_dat0$T2 - baseT)/12
# RMS_dat0$T3 <- (RMS_dat0$T3 - baseT)/12
# RMS_dat0$T4 <- (RMS_dat0$T4 - baseT)/12
# RMS_dat0$T5 <- (RMS_dat0$T5 - baseT)/12
# RMS_dat0$T6 <- (RMS_dat0$T6 - baseT)/12
# RMS_dat0$T7 <- (RMS_dat0$T7 - baseT)/12
# RMS_dat0$T8 <- (RMS_dat0$T8 - baseT)/12
# RMS_dat0$T9 <- (RMS_dat0$T9 - baseT)/12
# # Standardize time-invariant covariates (TICs)
# ## ex1 is standardized growth TIC in models
# RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
# # Standardize time-varying covariate (TVC)
# BL_mean <- mean(RMS_dat0[, "R1"])
# BL_var <- var(RMS_dat0[, "R1"])
# RMS_dat0$Rs1 <- (RMS_dat0$R1 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs2 <- (RMS_dat0$R2 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs3 <- (RMS_dat0$R3 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs4 <- (RMS_dat0$R4 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs5 <- (RMS_dat0$R5 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs6 <- (RMS_dat0$R6 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs7 <- (RMS_dat0$R7 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs8 <- (RMS_dat0$R8 - BL_mean)/sqrt(BL_var)
# RMS_dat0$Rs9 <- (RMS_dat0$R9 - BL_mean)/sqrt(BL_var)
# xstarts <- mean(baseT)

## ----message = FALSE, eval = FALSE--------------------------------------------
# set.seed(20191029)
# Math_TVC_BLS_f <- getTVCmodel(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE,
#   records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 0, growth_TIC = "ex1",
#   tries = 10
#   )
# set.seed(20191029)
# Math_TVCslp_BLS_f <- getTVCmodel(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE,
#   records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 1, growth_TIC = "ex1",
#   tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
getEstimateStats(est_in = Math_TVCslp_BLS_f@Estimates, CI_type = "Wald")
Figure1 <- getFigure(
  model = Math_TVC_BLS_f@mxOutput, sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
show(Figure1)
Figure2 <- getFigure(
  model = Math_TVCslp_BLS_f@mxOutput, sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
show(Figure2)

## ----message = FALSE, eval = FALSE--------------------------------------------
# set.seed(20191029)
# Math_TVCslp_BLS_r <- getTVCmodel(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#   records = 1:9, y_model = "LGCM", TVC = "R", decompose = 1, growth_TIC = "ex1",
#   tries = 10, paramOut = TRUE)
# set.seed(20191029)
# Math_TVCchg_BLS_r <- getTVCmodel(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#   records = 1:9, y_model = "LGCM", TVC = "R", decompose = 2, growth_TIC = "ex1",
#   tries = 10, paramOut = TRUE)
# set.seed(20191029)
# Math_TVCchgBL_BLS_r <- getTVCmodel(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#   records = 1:9, y_model = "LGCM", TVC = "R", decompose = 3, growth_TIC = "ex1",
#   tries = 10, paramOut = TRUE)

## -----------------------------------------------------------------------------
Figure3 <- getFigure(
  model = Math_TVCslp_BLS_r@mxOutput, sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
show(Figure3)
Figure4 <- getFigure(
  model = Math_TVCchg_BLS_r@mxOutput, sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
show(Figure4)
Figure5 <- getFigure(
  model = Math_TVCchgBL_BLS_r@mxOutput, sub_Model = "TVC", y_var = "M", curveFun = "BLS", 
  y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Mathematics"
)
show(Figure5)

