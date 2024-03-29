## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ---- message = FALSE---------------------------------------------------------
load(system.file("extdata", "getMGM_examples.RData", package = "nlpsem"))

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
#  xstarts <- mean(baseT)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  paraBLS_PLGCM.r <- c(
#    "Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_knot",
#    paste0("Y_psi", c("00", "01", "02", "11", "12", "22")), "Y_res",
#    "Z_mueta0", "Z_mueta1", "Z_mueta2", "Z_knot",
#    paste0("Z_psi", c("00", "01", "02", "11", "12", "22")), "Z_res",
#    paste0("YZ_psi", c("00", "10", "20", "01", "11", "21", "02", "12", "22")),
#    "YZ_res"
#    )
#  RM_PLGCM.r <- getMGM(
#    dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
#    intrinsic = FALSE, records = list(1:9, 1:9), y_model = "LGCM", res_scale = c(0.1, 0.1),
#    res_cor = 0.3, paramOut = TRUE, names = paraBLS_PLGCM.r
#    )

## -----------------------------------------------------------------------------
Figure1 <- getFigure(
  model = RM_PLGCM.r@mxOutput, sub_Model = "MGM", y_var = c("R", "M"), curveFun = "BLS", 
  y_model = "LGCM", t_var = c("T", "T"), records = list(1:9, 1:9), xstarts = xstarts, 
  xlab = "Month", outcome = c("Reading", "Mathematics")
)
show(Figure1)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  paraBLS_PLGCM_f <- c(
#    "Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_knot",
#    paste0("Y_psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")), "Y_res",
#    "Z_mueta0", "Z_mueta1", "Z_mueta2", "Z_knot",
#    paste0("Z_psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")), "Z_res",
#    paste0("YZ_psi", c(c("00", "10", "20", "g0", "01", "11", "21", "g1",
#                         "02", "12", "22", "g2", "0g", "1g", "2g", "gg"))),
#    "YZ_res"
#    )
#  RM_PLGCM.f <- getMGM(
#    dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
#    intrinsic = TRUE, records = list(1:9, 1:9), y_model = "LGCM", res_scale = c(0.1, 0.1),
#    res_cor = 0.3, paramOut = TRUE, names = paraBLS_PLGCM_f
#    )

## -----------------------------------------------------------------------------
Figure2 <- getFigure(
  model = RM_PLGCM.f@mxOutput, sub_Model = "MGM", y_var = c("R", "M"), curveFun = "BLS", 
  y_model = "LGCM", t_var = c("T", "T"), records = list(1:9, 1:9), xstarts = xstarts, 
  xlab = "Month", outcome = c("Reading", "Mathematics")
)
show(Figure2)

