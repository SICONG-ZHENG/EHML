#' Run glmboost Classification Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct glmboost
#'   classification models.
#' @param cutoff_method Method used to determine the classification cutoff.
#'   One of \code{"fixed"}, \code{"youden"}, or
#'   \code{"closest.topleft"}.
#' @param cutoff Fixed probability cutoff used when
#'   \code{cutoff_method = "fixed"}.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted glmboost
#'   models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunGlmboost(
#'   Object = Object,
#'   gene = gene,
#'   cutoff_method = "youden",
#'   cutoff = 0.5,
#'   seed = 123
#' )
#' }
RunGlmboost <- function(Object,
                        gene,
                        cutoff_method = c("fixed","youden","closest.topleft"),
                        cutoff=0.5,
                        seed=123){
  #######################################
  Glmboost_result <- list()
  Glmboost_result_F1 <- list()
  Glmboost_result_ACC <- list()
  Glmboost_result_Recall <- list()
  Glmboost_result_NPV <- list()
  Glmboost_result_PPV <- list()
  Glmboost_result_Sen <-list()
  Glmboost_result_Spe <- list()
  Glmboost_result_auc <- list()
  Glmboost_result_summary <- list()
  Glmboost_result_Youden <-list()
  Glmboost_result_kappa <- list()
  Glmboost_result_model<- list()
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
        "Model: Glmboost","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")

    train_exp <- train_data[,selection,drop=F]
    train_exp$labels <- factor(train_label$Label,levels = c(0,1))

    set.seed(seed)
    model <- mboost::glmboost(labels ~ .,data = train_exp,family = Binomial())
    cvm <- mboost::cvrisk(model, papply = lapply,
                  folds = cv(model.weights(model), type = "kfold"))
    model <- mboost::glmboost(labels ~ .,
                      data = train_exp,
                      family = Binomial(),
                      control = boost_control(mstop = max(mstop(cvm), 40)))
    # model
    train_result <- as.data.frame(stats::predict(model,newdata = train_exp, type = "response"))
    train_result$real_label <- factor(ifelse(train_exp$labels==1, "positive", "negative"))
    colnames(train_result) <- c("predict_p", "real_label")


    roc_obj <- pROC::roc(response = train_result$real_label,predictor = train_result$predict_p, levels = c("negative", "positive"),direction = "<")
    Youden_cut <-  pROC::coords(roc_obj,x = "best",best.method = "youden",ret = c("threshold", "sensitivity", "specificity"))
    AUC_cut <-  pROC::coords(roc_obj,x = "best",best.method = "closest.topleft",ret = c("threshold", "sensitivity", "specificity"))

    ###############cutoff:Youden_or_AUC############################################################################################################################
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
    name=paste0(selectionname," + Glmboost(cutoff:", round(cutoff_use,2), ")")
    Glmboost_result_F1[[name]] <- result_F1
    Glmboost_result_ACC[[name]] <- result_ACC
    Glmboost_result_Recall[[name]] <- result_Recall
    Glmboost_result_NPV[[name]] <- result_NPV
    Glmboost_result_PPV[[name]] <- result_PPV
    Glmboost_result_Sen[[name]] <- result_Sen
    Glmboost_result_Spe[[name]] <- result_Spe
    Glmboost_result_auc[[name]] <- result_auc
    Glmboost_result_Youden[[name]] <- result_Youden
    Glmboost_result_kappa[[name]] <- result_kappa
    Glmboost_result_model[[name]] <- model
    Glmboost_result_summary[[name]] <- all_result
  }
  Glmboost_result <- list(
    F1 = Glmboost_result_F1,
    ACC = Glmboost_result_ACC,
    RECALL = Glmboost_result_Recall,
    NPV = Glmboost_result_NPV,
    PPV = Glmboost_result_PPV,
    Sensitive=Glmboost_result_Sen,
    Specific=Glmboost_result_Spe,
    AUC=Glmboost_result_auc,
    Youden=Glmboost_result_Youden,
    Kappa=Glmboost_result_kappa,
    MODEL =Glmboost_result_model,
    SUMMARY =Glmboost_result_summary)
  return(Glmboost_result)
}
