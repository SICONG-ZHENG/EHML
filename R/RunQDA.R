#' Run Quadratic Discriminant Analysis
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct QDA models.
#'
#' @return A list containing model performance metrics, fitted QDA models,
#'   and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunQDA(
#'   Object = Object,
#'   gene = gene
#' )
#' }
RunQDA<- function(Object,
                  gene){
  #######################################
  QDA_result <- list()
  QDA_result_F1 <- list()
  QDA_result_ACC <- list()
  QDA_result_Recall <- list()
  QDA_result_NPV <- list()
  QDA_result_PPV <- list()
  QDA_result_Sen <-list()
  QDA_result_Spe <- list()
  QDA_result_auc <- list()
  QDA_result_summary <- list()
  QDA_result_Youden <-list()
  QDA_result_kappa <- list()
  QDA_result_model<- list()
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
        "Model: Quadratic discriminant analysis(QDA)","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")
    train_exp <- train_data[,selection,drop=F]
    train_exp$labels <- train_label$Label
    model <- MASS::qda(labels~., data = train_exp)
    # if (any(counts < p + 1)) stop("some group is too small for 'qda'"),
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
    name=paste0(selectionname," + QDA")
    QDA_result_F1[[name]] <- result_F1
    QDA_result_ACC[[name]] <- result_ACC
    QDA_result_Recall[[name]] <- result_Recall
    QDA_result_NPV[[name]] <- result_NPV
    QDA_result_PPV[[name]] <- result_PPV
    QDA_result_Sen[[name]] <- result_Sen
    QDA_result_Spe[[name]] <- result_Spe
    QDA_result_auc[[name]] <- result_auc
    QDA_result_Youden[[name]] <- result_Youden
    QDA_result_kappa[[name]] <- result_kappa
    QDA_result_model[[name]] <- model
    QDA_result_summary[[name]] <- all_result
  }
  QDA_result <- list( F1 = QDA_result_F1,
                      ACC = QDA_result_ACC,
                      RECALL = QDA_result_Recall,
                      NPV = QDA_result_NPV,
                      PPV = QDA_result_PPV,
                      Sensitive=QDA_result_Sen,
                      Specific=QDA_result_Spe,
                      AUC=QDA_result_auc,
                      Youden=QDA_result_Youden,
                      Kappa=QDA_result_kappa,
                      MODEL =QDA_result_model,
                      SUMMARY =QDA_result_summary)
  return(QDA_result)
}



