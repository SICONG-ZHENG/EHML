
#' Run k-Nearest Neighbor Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct KNN models.
#' @param knumber A numeric vector specifying the numbers of nearest
#'   neighbors to evaluate.
#'
#' @return A list containing model performance metrics, fitted KNN models,
#'   and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunKNN(
#'   Object = Object,
#'   gene = gene,
#'   knumber = c(1, 2, 3, 4, 5)
#' )
#' }
#'
#'
#'
RunKNN<- function(Object,
                  gene,
                  knumber = c(1, 2, 3, 4, 5)){
  #######################################
  KNN_result <- list()
  KNN_result_F1 <- list()
  KNN_result_ACC <- list()
  KNN_result_Recall <- list()
  KNN_result_NPV <- list()
  KNN_result_PPV <- list()
  KNN_result_Sen <-list()
  KNN_result_Spe <- list()
  KNN_result_auc <- list()
  KNN_result_summary <- list()
  KNN_result_Youden <-list()
  KNN_result_kappa <- list()
  KNN_result_model<- list()
  ####
  train_data <- Object$train_data
  train_label <- Object$train_label
  #######################################
  test_explist <-  Object$test_list
  test_label <- Object$test_label
  #########################################
  for (i in 1:length(gene)) {
    for (knumberi in knumber) {
      selection <- gene[[i]]
      selectionname <- names(gene)[i]
      cat("\n====================================\n",
          "Model: k-Nearest Neighbor(KNN)","\n",
          "Gene Set:", selectionname, "\n",
          "Gene Number:", length(selection), "\n",
          "====================================\n")
      train_exp <- train_data[,selection,drop=F]
      #train_exp$labels <- train_label$Label
      z <- class::knn(train_exp, train_exp, train_label$Label,
                      k = knumberi, prob = TRUE)
      train_result <- data.frame(
        predict_result = factor(ifelse(z == 1, "positive", "negative")),
        predict_p = ifelse(z == 1, attr(z, "prob"), 1 - attr(z, "prob")),
        real_label = factor(ifelse(train_label$Label == 1, "positive", "negative"))
      )
      rownames(train_result) <- rownames(train_exp)
      model <- list(train_data = train_exp,train_label = train_label$Label,k = knumberi)
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
        #zsc
        ztest <- class::knn(model$train_data,expdata,model$train_label,model$k,prob = TRUE)
        tresult <- data.frame(
          predict_result = factor(ifelse(ztest == 1, "positive", "negative")),
          predict_p = ifelse(ztest == 1,attr(ztest, "prob"),
                             1 - attr(ztest, "prob")),
          real_label = factor(ifelse(real_label_numeric == 1, "positive", "negative")),
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
      name=paste0(selectionname," + KNN(k:",knumberi, ")")
      KNN_result_F1[[name]] <- result_F1
      KNN_result_ACC[[name]] <- result_ACC
      KNN_result_Recall[[name]] <- result_Recall
      KNN_result_NPV[[name]] <- result_NPV
      KNN_result_PPV[[name]] <- result_PPV
      KNN_result_Sen[[name]] <- result_Sen
      KNN_result_Spe[[name]] <- result_Spe
      KNN_result_auc[[name]] <- result_auc
      KNN_result_Youden[[name]] <- result_Youden
      KNN_result_kappa[[name]] <- result_kappa
      KNN_result_model[[name]] <- model
      KNN_result_summary[[name]] <- all_result
    }}
  KNN_result <- list(F1 = KNN_result_F1,
                     ACC = KNN_result_ACC,
                     RECALL = KNN_result_Recall,
                     NPV = KNN_result_NPV,
                     PPV = KNN_result_PPV,
                     Sensitive=KNN_result_Sen,
                     Specific=KNN_result_Spe,
                     AUC=KNN_result_auc,
                     Youden=KNN_result_Youden,
                     Kappa=KNN_result_kappa,
                     MODEL =KNN_result_model,
                     SUMMARY =KNN_result_summary)
  return(KNN_result)
}

