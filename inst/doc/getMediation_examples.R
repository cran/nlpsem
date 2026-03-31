## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_data <- nzchar(system.file("extdata", "getMediation_examples.RData", package = "nlpsem"))
knitr::opts_chunk$set(eval = has_data)

## ----message = FALSE----------------------------------------------------------
library(nlpsem)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

## ----message = FALSE----------------------------------------------------------
load(system.file("extdata", "getMediation_examples.RData", package = "nlpsem"))

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
# RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
# xstarts <- mean(baseT)

## ----message = FALSE, eval = FALSE--------------------------------------------
# Med2_LGCM_BLS <- getMediation(
#   dat = RMS_dat0, t_var = rep("T", 2), y_var = "M", m_var = "R",
#   x_type = "baseline", x_var = "ex1", curveFun = "bilinear spline",
#   records = list(1:9, 1:9), tries = 10,
#   paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Med2_LGCM_BLS@Estimates

## ----message = FALSE, eval = FALSE--------------------------------------------
# set.seed(20191029)
# Med3_LGCM_BLS <- getMediation(
#   dat = RMS_dat0, t_var = rep("T", 3), y_var = "S", m_var = "M", x_type = "longitudinal",
#   x_var = "R", curveFun = "bilinear spline", records = list(2:9, 1:9, 1:9),
#   tries = 10, paramOut = TRUE
#   )

## -----------------------------------------------------------------------------
Med3_LGCM_BLS@Estimates

