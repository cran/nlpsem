#' @title Fit a Longitudinal Mediation Model
#'
#' @description This function fits a longitudinal mediation model to the provided data. It manages model setup, optimization,
#' and if requested, outputs parameter estimates and standard errors.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables
#' with repeated measurements and occasions for multiple longitudinal processes and a baseline predictor when applicable.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time
#' variable for the corresponding longitudinal variable at each study wave.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' @param m_var A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options include: \code{"linear"}
#' (or \code{"LIN"}), and \code{"bilinear spline"} (or \code{"BLS"}).
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding longitudinal variable.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified initial
#' values.
#' @param res_scale A numeric vector with each element representing the scaling factor for the initial calculation of the residual
#' variance. These values should be between \code{0} and \code{1}, exclusive. By default, this is \code{NULL}, as it is unnecessary
#' when the user specifies the initial values using the \code{starts} argument.
#' @param res_cor A numeric value or vector for user-specified residual correlation between any two longitudinal processes to calculate
#' the corresponding initial value. By default, this is \code{NULL}, as it is unnecessary when the user specifies the initial values
#' using the \code{starts} argument.
#' @param tries An integer specifying the number of additional optimization attempts. Default is \code{NULL}.
#' @param OKStatus An integer (vector) specifying acceptable status codes for convergence. Default is \code{0}.
#' @param jitterD A string specifying the distribution for jitter. Supported values are: \code{"runif"} (uniform
#' distribution), \code{"rnorm"} (normal distribution), and \code{"rcauchy"} (Cauchy distribution). Default is \code{"runif"}.
#' @param loc A numeric value representing the location parameter of the jitter distribution. Default is \code{1}.
#' @param scale A numeric value representing the scale parameter of the jitter distribution. Default is \code{0.25}.
#' @param paramOut A logical flag indicating whether to output the parameter estimates and standard errors. Default is \code{FALSE}.
#' @param names A character vector specifying parameter names. Default is \code{NULL}.
#'
#' @return An object of class \code{myMxOutput}. Depending on the \code{paramOut} argument, the object may contain the following slots:
#' \itemize{
#'   \item \code{mxOutput}: This slot contains the fitted longitudinal mediation model. A summary of this model can be obtained using
#'   the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R.A. (2022). Assessing Mediational Processes Using Piecewise Linear Growth Curve Models with Individual
#'   Measurement Occasions. Behavior Research Methods (Advance online publication). \doi{10.3758/s13428-022-01940-2}}
#'   \item {MacKinnon, D. P. (2008). Introduction to Statistical Mediation Analysis. Taylor & Francis Group/Lawrence Erlbaum Associates.}
#'   \item {Cheong, J., Mackinnon, D. P., & Khoo, S. T. (2003). Investigation of Mediational Processes Using Parallel Process Latent
#'   Growth Curve Modeling. Structural equation modeling: a multidisciplinary journal, 10(2), 238-262.
#'   \doi{10.1207/S15328007SEM1002_5}}
#'   \item {Soest, T., & Hagtvet, K. A. (2011). Mediation Analysis in a Latent Growth Curve Modeling Framework. Structural equation modeling:
#'   a multidisciplinary journal, 18(2), 289-314. \doi{10.1080/10705511.2011.557344}}
#' }
#'
#' @export
#'
#' @examples
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load ECLS-K (2011) data
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the starting point of the study
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#'
#' \donttest{
#' # Example 1: Baseline predictor, linear functional form
#' ## Fit model
#' set.seed(20191029)
#' Med2_LGCM_LIN <- getMediation(
#'   dat = RMS_dat0, t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline",
#'   x_var = "ex1", curveFun = "LIN", records = list(1:9, 1:9), res_scale = c(0.1, 0.1),
#'   res_cor = 0.3
#'   )
#'
#' # Example 2: Longitudinal predictor, bilinear spline functional form
#' ## Define parameter names
#' paraMed3_BLS <- c(
#'   "muetaX1", "muetaXr", "muetaX2", "mugX",
#'   paste0("psi", c("X1X1", "X1Xr", "X1X2", "XrXr", "XrX2", "X2X2")),
#'   "alphaM1", "alphaMr", "alphaM2", "mugM",
#'   paste0("psi", c("M1M1", "M1Mr", "M1M2", "MrMr", "MrM2", "M2M2"), "_r"),
#'   "alphaY1", "alphaYr", "alphaY2", "mugY",
#'   paste0("psi", c("Y1Y1", "Y1Yr", "Y1Y2", "YrYr", "YrY2", "Y2Y2"), "_r"),
#'   paste0("beta", c("X1Y1", "X1Yr", "X1Y2", "XrYr", "XrY2", "X2Y2",
#'                    "X1M1", "X1Mr", "X1M2", "XrMr", "XrM2", "X2M2",
#'                    "M1Y1", "M1Yr", "M1Y2", "MrYr", "MrY2", "M2Y2")),
#'   "muetaM1", "muetaMr", "muetaM2", "muetaY1", "muetaYr", "muetaY2",
#'   paste0("mediator", c("111", "11r", "112", "1rr", "1r2",
#'                        "122", "rr2", "r22", "rrr", "222")),
#'   paste0("total", c("11", "1r", "12", "rr", "r2", "22")),
#'   "residualsX", "residualsM", "residualsY", "residualsMX", "residualsYX", "residualsYM"
#'   )
#'
#' ## Fit model
#' set.seed(20191029)
#' Med3_LGCM_BLS <- getMediation(
#'   dat = RMS_dat0, t_var = rep("T", 3), y_var = "S", m_var = "M", x_type = "longitudinal",
#'   x_var = "R", curveFun = "bilinear spline", records = list(2:9, 1:9, 1:9),
#'   res_scale = c(0.1, 0.1, 0.1), res_cor = c(0.3, 0.3), tries = 10, paramOut = TRUE,
#'   names = paraMed3_BLS
#'   )
#' printTable(Med3_LGCM_BLS)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getMediation <- function(dat, t_var, y_var, m_var, x_type, x_var, curveFun, records, starts = NULL, res_scale = NULL,
                         res_cor = NULL, tries = NULL, OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25,
                         paramOut = FALSE, names = NULL){
  if (I(paramOut & is.null(names))){
    stop("Please enter the original parameters if want to obtain them!")
  }
  if (I(any(res_scale <= 0) | any(res_scale >= 1))){
    stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
  }
  if (!I(curveFun %in% c("linear", "LIN", "bilinear spline", "BLS"))){
    stop("Longitudinal mediation model only allows for linear or bilinear spline functional form!")
  }
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getMED.initial(dat = dat, t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type,
                             x_var = x_var, curveFun = curveFun, records = records, res_scale = res_scale,
                             res_cor = res_cor)
  }
  ## Build up a longitudinal mediation model
  ### Obtain manifest and latent variables, paths of the predictor, mediator and outcome
  model_mx <- getMED.mxModel(dat = dat, t_var = t_var, y_var = y_var,  m_var = m_var, x_type = x_type,
                             x_var = x_var, curveFun = curveFun, records = records, res_cor = res_cor,
                             starts = starts)
  if (!is.null(tries)){
    model0 <- mxTryHard(model_mx, extraTries = tries, OKstatuscodes = OKStatus, jitterDistrib = jitterD,
                       loc = loc, scale = scale)
    model <- mxRun(model0)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(paramOut){
    MED_output <- getMED.output(model = model, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                                curveFun = curveFun, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = MED_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  return(model)
}
