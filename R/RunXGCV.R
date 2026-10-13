#' Run XGBoost Models with Cross-Validation
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct XGBoost models.
#' @param nrounds Maximum number of boosting iterations evaluated during
#'   hyperparameter tuning.
#' @param fold Number of folds used for cross-validation.
#' @param cutoff_method Method used to determine the classification cutoff.
#'   One of \code{"fixed"}, \code{"youden"}, or
#'   \code{"closest.topleft"}.
#' @param cutoff Fixed probability cutoff used when
#'   \code{cutoff_method = "fixed"}.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted XGBoost
#'   models, and prediction results.
#' @export
#'
#' @details
#' This function relies on \pkg{caret} with the \code{"xgbTree"} method.
#' For compatibility with \pkg{caret}, xgboost version 1.7.8.1 is
#' recommended. Newer versions (>= 2.0.0) may introduce API changes that
#' are incompatible with the current implementation.
#'
#' @examples
#' \dontrun{
#' result <- RunXGCV(
#'   Object = Object,
#'   gene = gene,
#'   nrounds = 1000,
#'   fold = 10,
#'   cutoff_method = "youden",
#'   cutoff = 0.5,
#'   seed = 123
#' )
#' }
RunXGCV <- function(Object,
                    gene, nrounds= 1000,
                    fold=10,
                    cutoff_method = c("fixed","youden","closest.topleft"),
                    cutoff=0.5,
                    seed=123){
  #######################################
  XGCV_result <- list()
  XGCV_result_F1 <- list()
  XGCV_result_ACC <- list()
  XGCV_result_Recall <- list()
  XGCV_result_NPV <- list()
  XGCV_result_PPV <- list()
  XGCV_result_Sen <-list()
  XGCV_result_Spe <- list()
  XGCV_result_auc <- list()
  XGCV_result_summary <- list()
  XGCV_result_Youden <-list()
  XGCV_result_kappa <- list()
  XGCV_result_model<- list()
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
        "Model: XGBoost-CV","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")
    train_exp <- as.data.frame(train_data[,selection,drop=F])
    train_exp$labels <- factor(train_label$Label)

    set.seed(seed)
    tune_grid <- expand.grid(nrounds = seq(from = 50, to = nrounds, by = 50),
                             eta = c(0.025, 0.05, 0.1, 0.3),
                             max_depth = c(2, 3, 4, 5, 6),
                             gamma = 0,
                             colsample_bytree = 1,
                             min_child_weight = 1,
                             subsample = 1)
    fitControl <- caret::trainControl(method = "cv",
                                      number = fold,
                                      verboseIter = FALSE,
                                      search = "grid")
    set.seed(seed)
    caret_xgb <- caret::train(labels~. ,
                              data = train_exp,
                              method="xgbTree",
                              trControl=fitControl,
                              tuneGrid = tune_grid)
    model <- caret_xgb$finalModel
    # 评估模型
    train_result <- as.data.frame(predict(model, as.matrix(train_data[,selection,drop=F])))
    train_result$real_label <- factor(ifelse(train_label$Label==1, "positive", "negative"))
    names(train_result) <- c("predict_p","real_label")
    roc_obj <- pROC::roc(response = train_result$real_label,predictor = train_result$predict_p, levels = c("negative", "positive"),direction = "<")
    Youden_cut <-  pROC::coords(roc_obj,x = "best",best.method = "youden",ret = c("threshold", "sensitivity", "specificity"))
    AUC_cut <-  pROC::coords(roc_obj,x = "best",best.method = "closest.topleft",ret = c("threshold", "sensitivity", "specificity"))
    ###############cutoff:#Youden#or#AUC############################################################################################################################
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
    train_result$type <- factor(ifelse(train_result$predict_p > cutoff_use, "positive", "negative") )
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
      predict_p <-predict(model, as.matrix(expdata))
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
    name=paste0(selectionname," + XGBoost-CV(cutoff:", round(cutoff_use,2), ")")
    XGCV_result_F1[[name]] <- result_F1
    XGCV_result_ACC[[name]] <- result_ACC
    XGCV_result_Recall[[name]] <- result_Recall
    XGCV_result_NPV[[name]] <- result_NPV
    XGCV_result_PPV[[name]] <- result_PPV
    XGCV_result_Sen[[name]] <- result_Sen
    XGCV_result_Spe[[name]] <- result_Spe
    XGCV_result_auc[[name]] <- result_auc
    XGCV_result_Youden[[name]] <- result_Youden
    XGCV_result_kappa[[name]] <- result_kappa
    XGCV_result_model[[name]] <- model
    XGCV_result_summary[[name]] <- all_result
  }
  XGCV_result <- list( F1 = XGCV_result_F1,
                       ACC = XGCV_result_ACC,
                       RECALL = XGCV_result_Recall,
                       NPV = XGCV_result_NPV,
                       PPV = XGCV_result_PPV,
                       Sensitive=XGCV_result_Sen,
                       Specific=XGCV_result_Spe,
                       AUC=XGCV_result_auc,
                       Youden=XGCV_result_Youden,
                       Kappa=XGCV_result_kappa,
                       MODEL =XGCV_result_model,
                       SUMMARY =XGCV_result_summary)
  return(XGCV_result)
}
