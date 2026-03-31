## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getMGM_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getMGM_examples.RData", package = "nlpsem"))

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
# xstarts <- mean(baseT)

## ----message = FALSE, eval = FALSE--------------------------------------------
# RM_PLGCM.r <- getMGM(
#   dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
#   intrinsic = FALSE, records = list(1:9, 1:9), y_model = "LGCM",
#   tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Figure1 <- getFigure(
  model = RM_PLGCM.r@mxOutput, sub_Model = "MGM", y_var = c("R", "M"), curveFun = "BLS", 
  y_model = "LGCM", t_var = c("T", "T"), records = list(1:9, 1:9), xstarts = xstarts, 
  xlab = "Month", outcome = c("Reading", "Mathematics")
)
show(Figure1)

## ----message = FALSE, eval = FALSE--------------------------------------------
# RM_PLGCM.f <- getMGM(
#   dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
#   intrinsic = TRUE, records = list(1:9, 1:9), y_model = "LGCM",
#   tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Figure2 <- getFigure(
  model = RM_PLGCM.f@mxOutput, sub_Model = "MGM", y_var = c("R", "M"), curveFun = "BLS", 
  y_model = "LGCM", t_var = c("T", "T"), records = list(1:9, 1:9), xstarts = xstarts, 
  xlab = "Month", outcome = c("Reading", "Mathematics")
)
show(Figure2)

