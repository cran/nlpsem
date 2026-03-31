## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# library(nlpsem)
# mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
# 
# # Load ECLS-K (2011) data
# data("RMS_dat")
# RMS_dat0 <- RMS_dat
# 
# # Re-baseline the data so that the estimated initial status corresponds to
# # the starting point of the study
# baseT <- RMS_dat0$T1
# for (i in 1:9) {
#   RMS_dat0[[paste0("T", i)]] <- RMS_dat0[[paste0("T", i)]] - baseT
# }
# xstarts <- mean(baseT)

## ----eval = FALSE-------------------------------------------------------------
# # Fit bilinear spline LGCM (fixed knot)
# Math_BLS_fixed <- getLGCM(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#   intrinsic = FALSE, records = 1:9
# )

## ----eval = FALSE-------------------------------------------------------------
# # Fit with parameter output
# Math_BLS_fixed <- getLGCM(
#   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#   intrinsic = FALSE, records = 1:9,
#   paramOut = TRUE
# )
# 
# # View estimates
# printTable(Math_BLS_fixed)

## ----eval = FALSE-------------------------------------------------------------
# getSummary(model_list = list(Math_BLS_fixed@mxOutput))

## ----eval = FALSE-------------------------------------------------------------
# getEstimateStats(
#   est_in = Math_BLS_fixed@Estimates, CI_type = "Wald"
# )

## ----eval = FALSE-------------------------------------------------------------
# Figure <- getFigure(
#   model = Math_BLS_fixed@mxOutput, nClass = NULL, cluster_TIC = NULL,
#   sub_Model = "LGCM", y_var = "M", curveFun = "BLS", y_model = "LGCM",
#   t_var = "T", records = 1:9, m_var = NULL, x_var = NULL, x_type = NULL,
#   xstarts = xstarts, xlab = "Month", outcome = "Mathematics"
# )
# show(Figure)

