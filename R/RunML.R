#' Run Multiple Machine Learning Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct machine learning models.
#' @param method A character vector specifying the machine learning methods to run.
#'   Available methods include \code{"LR"}, \code{"LDA"}, \code{"QDA"},
#'   \code{"KNN"}, \code{"DT"}, \code{"RF"}, \code{"XG"}, \code{"XGCV"},
#'   \code{"RRCV"}, \code{"LACV"}, \code{"EnetCV"}, \code{"SVM"},
#'   \code{"SVMCV"}, \code{"GBM"}, \code{"GBMCV"}, \code{"NV"},
#'   \code{"RSF"}, \code{"Glmboost"}, \code{"Plsrglm"}, and \code{"NN"}.
#' @param seed Random seed used by models that support randomization.
#' @param fold Number of folds used by models that support cross-validation.
#' @param cutoff_method Common method used to determine the classification cutoff
#'   for models that support cutoff selection. Available options include
#'   \code{"fixed"}, \code{"youden"}, and \code{"closest.topleft"}.
#' @param cutoff Common fixed classification cutoff. If \code{NULL}, the default
#'   cutoff defined in each model-specific function is used.
#' @param control A named list of model-specific parameter settings. Parameters
#'   supplied in \code{control} override the corresponding common parameters.
#'   Parameters not specified use the defaults defined in the corresponding
#'   model-specific function.
#'
#' @return A merged list containing performance metrics, fitted models,
#'   and prediction results from all selected machine learning methods.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunML(
#'   Object = Object,
#'   gene = gene,
#'   method = c("LR", "LDA", "RF", "LACV"),
#'   seed = 123,
#'   fold = 10,
#'   cutoff_method = "youden"
#' )
#'
#' result <- RunML(
#'   Object = Object,
#'   gene = gene,
#'   method = c("KNN", "LACV", "NN"),
#'   control = list(
#'     KNN = list(
#'       knumber = c(3, 5, 7)
#'     ),
#'     LACV = list(
#'       fold = 5,
#'       cutoff_method = "fixed",
#'       cutoff = 0.5
#'     ),
#'     NN = list(
#'       learningrate = c(0.001, 0.01),
#'       dropoutrate = c(0.25, 0.5),
#'       epochselecti = 50
#'     )
#'   )
#' )
#' }
RunML <- function(
    Object,
    gene,
    method = c("LR","LDA","QDA","KNN","DT","RF","XG","XGCV",
               "RRCV","LACV","EnetCV","SVM","SVMCV","GBM",
               "GBMCV","NV","RSF","Glmboost","Plsrglm","NN"),
    seed = 123,
    fold = 10,
    cutoff_method = "youden",
    cutoff = NULL,
    control = list()
) {
  ############################################################
  # Method -> function
  ############################################################
  method_fun <- c(LR       = "RunLR",LDA      = "RunLDA",QDA      = "RunQDA",
                  KNN      = "RunKNN",DT       = "RunDT",RF       = "RunRF",
                  XG       = "RunXG",XGCV     = "RunXGCV",RRCV     = "RunRRCV",
                  LACV     = "RunLACV",EnetCV   = "RunEnetCV",SVM      = "RunSVM",
                  SVMCV    = "RunSVMCV",GBM      = "RunGBM",GBMCV    = "RunGBMCV",
                  NV       = "RunNV",RSF      = "RunRSF",Glmboost = "RunGlmboost",
                  Plsrglm  = "RunPlsrglm",NN       = "RunNN")
  ############################################################
  # Check method
  ############################################################
  wrong_method <- setdiff(
    method,
    names(method_fun))
  if (length(wrong_method) > 0) {
    stop("Unknown method: ",paste(wrong_method, collapse = ", "))
  }
  ############################################################
  # Common parameters
  ############################################################
  common_control <- list(seed = seed,fold = fold,
                         cutoff_method = cutoff_method,cutoff = cutoff)
  # NULL means: use RunXXX default value
  common_control <- common_control[
    !vapply(common_control,is.null,logical(1)
    )]
  ############################################################
  # Run
  ############################################################
  result_list <- list()
  for (m in method) {
    cat(
      "\n====================================\n",
      "Running:", m, "\n",
      "====================================\n")
    ##########################################################
    # Get function
    ##########################################################
    fun_name <- method_fun[[m]]
    if (!exists(fun_name, mode = "function")) {
      stop("Function '",fun_name,"' does not exist.")
    }
    fun <- get(fun_name,mode = "function")
    valid_args <- names(formals(fun))
    ##########################################################
    # Basic arguments
    ##########################################################
    args <- list(Object = Object,gene = gene)
    ##########################################################
    # Add common parameters only when supported
    ##########################################################
    common_use <- common_control[names(common_control) %in% valid_args]
    args <- c(args,common_use)
    ##########################################################
    # Model-specific control
    ##########################################################
    method_control <- control[[m]]
    if (!is.null(method_control)) {
      if (!is.list(method_control)) {
        stop("control[['",m,"']] must be a list.")}
      wrong_args <- setdiff(names(method_control), valid_args)
      if (length(wrong_args) > 0) {
        warning(m," does not contain parameter(s): ",
                paste(wrong_args, collapse = ", "),
                ". These parameters were ignored." )
        method_control <- method_control[
          names(method_control) %in% valid_args]
      }
      ########################################################
      # Model-specific parameter overrides common parameter
      ########################################################
      args[names(method_control)] <- method_control
    }
    ##########################################################
    # Run model
    ##########################################################
    result_list[[m]] <- do.call(fun,args
    )
  }
  ############################################################
  # Merge
  ############################################################
  result <- do.call(MergeML,result_list)
  return(result)
}
