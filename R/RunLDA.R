#' Run Linear Discriminant Analysis
#'
#' Fit LDA models using different feature sets and evaluate performance
#' in the training and external validation datasets.
#'
#' @param Object A list containing:
#'   \code{train_data}, \code{train_label},
#'   \code{test_list}, and \code{test_label}.
#' @param gene A named list of feature vectors.
#'
#' @return A list containing performance metrics, fitted models,
#'   and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunLDA(
#'   Object = data_object,
#'   gene = gene_list
#' )
#' }
RunLDA<- function(Object,
                  gene){
  #######################################
  LDA_result <- list()
  LDA_result_F1 <- list()
  LDA_result_ACC <- list()
  LDA_result_Recall <- list()
  LDA_result_NPV <- list()
  LDA_result_PPV <- list()
  LDA_result_Sen <-list()
  LDA_result_Spe <- list()
  LDA_result_auc <- list()
  LDA_result_summary <- list()
  LDA_result_Youden <-list()
  LDA_result_kappa <- list()
  LDA_result_model<- list()
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
        "Model: Linear discriminant analysis(LDA)","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")
    train_exp <- train_data[,selection,drop=F]
    train_exp$labels <- train_label$Label
    model <- MASS::lda(labels~., data = train_exp)
    #zsc
    pretrain <- predict(model, type="response")
    train_result <- data.frame(predict_class=pretrain$class,
                               predict_result=pretrain$posterior[,2])
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("predict_class", "predict_p", "predict_result")
    train_result$real_label <- factor(ifelse(train_label==1, "positive", "negative"))
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
      predict_p <- predict(model, expdata, type="response")
      tresult <- data.frame(
        predict_p = predict_p$posterior[,2],
        real_label = factor(
          ifelse(real_label_numeric == 1,"positive","negative"),
          levels = c("negative", "positive")),
        predict_result = factor(
          ifelse(predict_p $class==1,"positive","negative"),
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
    name=paste0(selectionname," + LDA")
    LDA_result_F1[[name]] <- result_F1
    LDA_result_ACC[[name]] <- result_ACC
    LDA_result_Recall[[name]] <- result_Recall
    LDA_result_NPV[[name]] <- result_NPV
    LDA_result_PPV[[name]] <- result_PPV
    LDA_result_Sen[[name]] <- result_Sen
    LDA_result_Spe[[name]] <- result_Spe
    LDA_result_auc[[name]] <- result_auc
    LDA_result_Youden[[name]] <- result_Youden
    LDA_result_kappa[[name]] <- result_kappa
    LDA_result_model[[name]] <- model
    LDA_result_summary[[name]] <- all_result

  }
  LDA_result <- list( F1 = LDA_result_F1,
                      ACC = LDA_result_ACC,
                      RECALL = LDA_result_Recall,
                      NPV = LDA_result_NPV,
                      PPV = LDA_result_PPV,
                      Sensitive=LDA_result_Sen,
                      Specific=LDA_result_Spe,
                      AUC=LDA_result_auc,
                      Youden=LDA_result_Youden,
                      Kappa=LDA_result_kappa,
                      MODEL =LDA_result_model,
                      SUMMARY =LDA_result_summary)
  return(LDA_result)
}


