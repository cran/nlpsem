## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getLCSM_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getLCSM_examples.RData", package = "nlpsem"))

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
# ## ex1 and ex2 are standardized growth TICs in models
# RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
# RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
# xstarts <- mean(baseT)/12

## ----message = FALSE, eval = FALSE--------------------------------------------
# Read_LCSM_NonP <- getLCSM(
#   dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "nonparametric",
#   intrinsic = FALSE, records = 1:9, growth_TIC = NULL,
#   paramOut = TRUE
#   )
# Read_LCSM_NonP_TIC <- getLCSM(
#   dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "nonparametric",
#   intrinsic = FALSE, records = 1:9, growth_TIC = c("ex1", "ex2"),
#   paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
getSummary(model_list = list(Read_LCSM_NonP@mxOutput, Read_LCSM_NonP_TIC@mxOutput))
Figure1 <- getFigure(
  model = Read_LCSM_NonP@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "NonP", 
  y_model = "LCSM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Reading"
)
show(Figure1)
Figure2 <- getFigure(
  model = Read_LCSM_NonP_TIC@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "NonP", 
  y_model = "LCSM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Reading"
)
show(Figure2)

## ----message = FALSE, eval = FALSE--------------------------------------------
# Read_LCSM_QUAD <- getLCSM(
#   dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "quadratic", intrinsic = FALSE,
#   records = 1:9, paramOut = TRUE
#   )
# set.seed(20191029)
# Read_LCSM_EXP_r <- getLCSM(
#   dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "negative exponential",
#   intrinsic = FALSE, records = 1:9, tries = 10, paramOut = TRUE
#   )
# set.seed(20191029)
# Read_LCSM_JB_r <- getLCSM(
#   dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "Jenss-Bayley",
#   intrinsic = FALSE, records = 1:9, tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Figure3 <- getFigure(
  model = Read_LCSM_QUAD@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "QUAD", 
  y_model = "LCSM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Reading"
)
show(Figure3)
Figure4 <- getFigure(
  model = Read_LCSM_EXP_r@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "EXP", 
  y_model = "LCSM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Reading"
)
show(Figure4)
Figure5 <- getFigure(
  model = Read_LCSM_JB_r@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "JB", 
  y_model = "LCSM", t_var = "T", records = 1:9, xstarts = xstarts, xlab = "Year",
  outcome = "Reading"
)
show(Figure5)

