
#'  Run Logistic Regression Models
#'@description
#' Fits and evaluates binary logistic regression models using a training
#' dataset and multiple validation datasets. Models are constructed
#' separately for each supplied feature set.
#'
#' @param Object A list containing `train_data`, `train_label`,
#'   `test_list`, and `test_label`. The training and validation labels
#'   must contain a column named `Label` with binary values 0 and 1.
#' @param gene A named list of character vectors. Each element contains
#'   the predictor names used to construct a logistic regression model.
#' @param fold A positive integer indicating the number of
#'   cross-validation folds. Currently reserved and not used by
#'   this function.
#' @param cutoff_method Character string specifying the classification
#'   threshold selection method. Available options are `"fixed"`,
#'   `"youden"`, and `"closest.topleft"`.
#' @param cutoff Numeric classification threshold between 0 and 1,
#'   used when `cutoff_method = "fixed"`. Defaults to 0.5.
#' @param seed Integer random seed. Currently reserved and not used
#'   by this function.
#'
#' @return A named list containing the following components:
#'   \describe{
#'     \item{F1}{F1 scores for each model and dataset.}
#'     \item{ACC}{Classification accuracy.}
#'     \item{RECALL}{Recall scores.}
#'     \item{NPV}{Negative predictive values.}
#'     \item{PPV}{Positive predictive values.}
#'     \item{Sensitive}{Sensitivity.}
#'     \item{Specific}{Specificity.}
#'     \item{AUC}{Area under the ROC curve.}
#'     \item{Youden}{Youden indices.}
#'     \item{Kappa}{Cohen's kappa coefficients.}
#'     \item{MODEL}{Fitted logistic regression models.}
#'     \item{SUMMARY}{Sample-level predicted probabilities,
#'       true labels, and predicted classifications.}
#'   }
#'
#' @details
#' Logistic regression models are fitted using `stats::glm()`
#' with a binomial family.
#'
#' For each feature set, the classification threshold is determined
#' using either a fixed value, the Youden index, or the ROC point
#' closest to the top-left corner.
#'
#' Threshold optimization is performed using the training dataset.
#' The selected threshold is subsequently applied to the validation
#' datasets.
#'
#' Model performance is evaluated using multiple classification
#' metrics, including AUC, accuracy, F1 score, sensitivity,
#' specificity, predictive values, Youden index, and Cohen's kappa.
#'
#' @export
#'
RunLR <- function(Object,
                  gene,
                  fold=10,
                  cutoff_method = c("fixed","youden","closest.topleft"),
                  cutoff=0.5,
                  seed=123){
  #######################################
  LR_result <- list()
  LR_result_F1 <- list()
  LR_result_ACC <- list()
  LR_result_Recall <- list()
  LR_result_NPV <- list()
  LR_result_PPV <- list()
  LR_result_Sen <-list()
  LR_result_Spe <- list()
  LR_result_auc <- list()
  LR_result_summary <- list()
  LR_result_Youden <-list()
  LR_result_kappa <- list()
  LR_result_model<- list()
  #########
  cutoff_method <- match.arg(cutoff_method)
  if (length(cutoff) != 1 ||is.na(cutoff) ||!is.numeric(cutoff) ||!is.finite(cutoff) ||cutoff < 0 ||cutoff > 1) {
    warning("Invalid cutoff value. The default value (0.5) will be used.")
    cutoff <- 0.5}
  ####
  train_data <- Object$train_data
  train_label <- Object$train_label
  #######################################
  test_explist <-  Object$test_list
  test_label <- Object$test_label
#########################################
  for (i in 1:length(gene)) {
    selection <- gene[[i]]
    selectionname <- names(gene)[i]
    cat("\n====================================\n",
        "Model: Logistic Regression(LR)","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")
    train_exp <- train_data[,selection,drop=F]
    train_exp$labels <- train_label$Label
    model <- glm(labels~.-labels, family = "binomial", data = train_exp)
    # model
    train_result <- as.data.frame(predict(model, type="response"))
    colnames(train_result) <- "predict_p"
    train_result$real_label <- factor(ifelse(train_exp$labels==1, "positive", "negative"))
    roc_obj <- pROC::roc(response = train_result$real_label,predictor = train_result$predict_p, levels = c("negative", "positive"),direction = "<")
    Youden_cut <-  pROC::coords(roc_obj,x = "best",best.method = "youden",ret = c("threshold", "sensitivity", "specificity"))
    AUC_cut <-  pROC::coords(roc_obj,x = "best",best.method = "closest.topleft",ret = c("threshold", "sensitivity", "specificity"))
 ###############cutoff:#Youden#or#AUC#########################################################################################################################
    cutoff1 <- as.numeric(Youden_cut$threshold)
    cutoff2 <- as.numeric(AUC_cut$threshold)
    #######################################################################################################################################################
    cutoff_use <- switch(cutoff_method,
                         fixed = cutoff,
                         youden = cutoff1,
                         closest.topleft = cutoff2)
    if (length(cutoff_use) != 1 ||is.na(cutoff_use) ||!is.finite(cutoff_use)) {
      warning("Gene set '",selectionname,"' cannot obtain a valid cutoff; ","the default value 0.5 is used.")
      cutoff_use <- 0.5
    }
    train_result$type <- factor(ifelse(train_result[,1] > cutoff_use, "positive", "negative") )
    colnames(train_result) <- c("predict_p", "real_label","predict_result")
    all_result <- Map(function(data, labelsdata) {

      data <- as.data.frame(data)
      labelsdata <- as.data.frame(labelsdata)
      comd <- intersect(rownames(data),rownames(labelsdata))
      if (length(comd) == 0) {
        warning("test set no paried sample")
        return(NULL)
      }
      expdata <- data[comd,selection,drop = FALSE]
      real_label_numeric <- labelsdata[comd,"Label"]
      predict_p <- as.numeric(stats::predict(model,newdata = expdata,type = "response"))
      tresult <- data.frame(
        predict_p = predict_p,
        real_label = factor(
          ifelse(real_label_numeric == 1,"positive","negative"),
          levels = c("negative", "positive")),
        predict_result = factor(
          ifelse(predict_p > cutoff_use,
                 "positive",
                 "negative"),
          levels = c("negative", "positive")),
        row.names = comd)
      return(tresult)},data = test_explist,labelsdata = test_label)
    all_result[["Train"]] <- train_result
    result_Youden <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      Sensitivity <- cm$byClass["Sensitivity"]
      Specificity <- cm$byClass["Specificity"]
      Youden =Sensitivity + Specificity - 1
      Youden =unname(Youden)
      return(Youden)})
    result_kappa <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      kappa =unname(cm$overall["Kappa"])
      return(kappa)})
    result_F1 <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      F1 =unname(cm$byClass["F1"])
      return(F1)})
    result_ACC <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      ACC=unname(cm$overall["Accuracy"])
      return(ACC)})
    result_Recall <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(
        data = x$predict_result,
        reference = x$real_label,
        positive = "positive",
        mode = "everything")
      Recall=unname(cm$byClass["Recall"])
      return(Recall)})
    result_NPV <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      NPV = unname(cm$byClass["Neg Pred Value"])
      return(NPV)})
    result_PPV <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      PPV = unname(cm$byClass["Pos Pred Value"])
      return(PPV)})
    result_Sen <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      Sensitivity = unname(cm$byClass["Sensitivity"])
      return(Sensitivity)})
    result_Spe <-sapply(all_result, function(x){
      x$predict_result <- factor(x$predict_result,levels = c("negative","positive"))
      x$real_label <- factor(x$real_label,levels = c("negative","positive"))
      cm <- caret::confusionMatrix(data = x$predict_result,
                                   reference = x$real_label,
                                   positive = "positive",
                                   mode = "everything")
      Specificity = unname(cm$byClass["Specificity"])
      return(Specificity)})
    result_auc <- sapply(all_result, function(x) {
      pROC::roc(response = ifelse(x$real_label == "positive", 1, 0),
                predictor = x$predict_p,
                quiet = TRUE)$auc})
    name=paste0(selectionname," + LR(cutoff:",round(cutoff_use, 2), ")")
    LR_result_F1[[name]] <- result_F1
    LR_result_ACC[[name]] <- result_ACC
    LR_result_Recall[[name]] <- result_Recall
    LR_result_NPV[[name]] <- result_NPV
    LR_result_PPV[[name]] <- result_PPV
    LR_result_Sen[[name]] <- result_Sen
    LR_result_Spe[[name]] <- result_Spe
    LR_result_auc[[name]] <- result_auc
    LR_result_Youden[[name]] <- result_Youden
    LR_result_kappa[[name]] <- result_kappa
    LR_result_model[[name]] <- model
    LR_result_summary[[name]] <- all_result

  }
  LR_result <- list(
    F1 = LR_result_F1,
    ACC = LR_result_ACC,
    RECALL = LR_result_Recall,
    NPV = LR_result_NPV,
    PPV = LR_result_PPV,
    Sensitive=LR_result_Sen,
    Specific=LR_result_Spe,
    AUC=LR_result_auc,
    Youden=LR_result_Youden,
    Kappa=LR_result_kappa,
    MODEL =LR_result_model,
    SUMMARY =LR_result_summary)
  return(LR_result)
}


