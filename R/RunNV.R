#' Run Naive Bayes Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct naive Bayes models.
#'
#' @return A list containing model performance metrics, fitted naive Bayes
#'   models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunNV(
#'   Object = Object,
#'   gene = gene
#' )
#' }
RunNV <- function(Object,
                  gene){
  #######################################
  NV_result <- list()
  NV_result_F1 <- list()
  NV_result_ACC <- list()
  NV_result_Recall <- list()
  NV_result_NPV <- list()
  NV_result_PPV <- list()
  NV_result_Sen <-list()
  NV_result_Spe <- list()
  NV_result_auc <- list()
  NV_result_summary <- list()
  NV_result_Youden <-list()
  NV_result_kappa <- list()
  NV_result_model<- list()
  #########

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
        "Model: NaiveBayes","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")
    train_exp <- train_data[,selection,drop=F]
    labels <- factor(train_label$Label,levels = c(0,1))
    model <- e1071::naiveBayes(x = train_exp,y = labels)
    train_result <- as.data.frame(predict(model, train_exp))
    rownames(train_result) <- rownames(train_exp)
    prob <- predict(model, train_exp, type = "raw")

    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    train_result$predict_p <- prob[,2]
    train_result$real_label <- factor(ifelse(labels==1, "positive", "negative"))
    train_result <- train_result[,c("predict_p","real_label","type")]
    colnames(train_result) <- c("predict_p", "real_label","predict_result")
    all_result <- Map(function(data, labelsdata) {
      data = as.data.frame(data)
      labelsdata <- as.data.frame(labelsdata)
      comd <- intersect(rownames(data),rownames(labelsdata))
      if (length(comd) == 0) {
        warning("test set no paried sample")
        return(NULL)
      }
      expdata <- data[comd,selection,drop = FALSE]
      real_label_numeric <- labelsdata[comd,"Label"]
      tpred <- as.data.frame(predict(model, expdata))
      #zsc
      probt <- predict(model, expdata, type = "raw")
      tresult <- data.frame(
        predict_p = probt[,2],
        real_label = factor(
          ifelse(real_label_numeric == 1,"positive","negative"),
          levels = c("negative", "positive")),
        predict_result = factor(
          ifelse(tpred[,1] == 1,"positive","negative"),
          levels = c("negative", "positive")),
        row.names = comd)
      return(tresult)
    },data = test_explist,labelsdata = test_label)
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
    name=paste0(selectionname," + NaiveBayes")
    NV_result_F1[[name]] <- result_F1
    NV_result_ACC[[name]] <- result_ACC
    NV_result_Recall[[name]] <- result_Recall
    NV_result_NPV[[name]] <- result_NPV
    NV_result_PPV[[name]] <- result_PPV
    NV_result_Sen[[name]] <- result_Sen
    NV_result_Spe[[name]] <- result_Spe
    NV_result_auc[[name]] <- result_auc
    NV_result_Youden[[name]] <- result_Youden
    NV_result_kappa[[name]] <- result_kappa
    NV_result_model[[name]] <- model
    NV_result_summary[[name]] <- all_result
  }
  NV_result <- list( F1 = NV_result_F1,
                     ACC = NV_result_ACC,
                     RECALL = NV_result_Recall,
                     NPV = NV_result_NPV,
                     PPV = NV_result_PPV,
                     Sensitive=NV_result_Sen,
                     Specific=NV_result_Spe,
                     AUC=NV_result_auc,
                     Youden=NV_result_Youden,
                     Kappa=NV_result_kappa,
                     MODEL =NV_result_model,
                     SUMMARY =NV_result_summary)
  return(NV_result)
}
