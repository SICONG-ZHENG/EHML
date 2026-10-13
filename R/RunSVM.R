#' Run Support Vector Machine Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct support vector
#'   machine models.
#' @param kernel_all A character vector specifying the kernel functions to
#'   evaluate. Available options include \code{"linear"},
#'   \code{"polynomial"}, and \code{"radial"}.
#' @param fold Number of cross-validation folds. Currently not used in
#'   the function.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted support
#'   vector machine models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunSVM(
#'   Object = Object,
#'   gene = gene,
#'   kernel_all = c("linear", "polynomial", "radial"),
#'   seed = 123
#' )
#' }
RunSVM<- function(Object,
                  gene,
                  kernel_all = c("linear", "polynomial", "radial"),
                  fold=10,
                  seed=123){
  #######################################
  SVM_result <- list()
  SVM_result_F1 <- list()
  SVM_result_ACC <- list()
  SVM_result_Recall <- list()
  SVM_result_NPV <- list()
  SVM_result_PPV <- list()
  SVM_result_Sen <-list()
  SVM_result_Spe <- list()
  SVM_result_auc <- list()
  SVM_result_summary <- list()
  SVM_result_Youden <-list()
  SVM_result_kappa <- list()
  SVM_result_model<- list()
  ####
  train_data <- Object$train_data
  train_label <- Object$train_label
  #######################################
  test_explist <-  Object$test_list
  test_label <- Object$test_label
  #########################################
  for (i in 1:length(gene)) {
    for (kerneli in kernel_all) {
      selection <- gene[[i]]
      selectionname <- names(gene)[i]
      cat("\n====================================\n",
          "Model: Support vector machine","\n",
          "Gene Set:", selectionname, "\n",
          "Gene Number:", length(selection), "\n",
          "Kernel:",kerneli, "\n",
          "====================================\n")
      train_exp <- train_data[,selection,drop=F]
      train_exp$labels <- train_label$Label
      train_exp$labels <- factor(train_exp$labels,levels = c(0,1))
      set.seed(seed)
      model<-e1071::svm(labels~. , data = train_exp,
                        kernel=kerneli,probability = TRUE)
      pred <- predict(model,as.data.frame(train_exp),probability = TRUE)
      prob <- attr(pred, "probabilities")
      train_result <- data.frame(real_label=factor(ifelse(train_exp$labels==1, "positive", "negative")),
                                 predict_result=predict(model,as.data.frame(train_exp),probability = TRUE),
                                 predict_p = prob[,colnames(prob)%in%1])
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
        predt <- predict(model,as.data.frame(expdata),probability = TRUE)
        probt <- attr(predt, "probabilities")
        tresult <- as.data.frame(predict(model, as.data.frame(expdata)))
        tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
        colnames(tresult) <- c("real_label", "predict_result")
        tresult$real_label <- factor(ifelse(real_label_numeric==1, "positive", "negative"))
        #zsc
        tresult$predict_p <- probt[,2]
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
      name=paste0(selectionname," + SVM-default(kernel:",kerneli,")")
      SVM_result_F1[[name]] <- result_F1
      SVM_result_ACC[[name]] <- result_ACC
      SVM_result_Recall[[name]] <- result_Recall
      SVM_result_NPV[[name]] <- result_NPV
      SVM_result_PPV[[name]] <- result_PPV
      SVM_result_Sen[[name]] <- result_Sen
      SVM_result_Spe[[name]] <- result_Spe
      SVM_result_auc[[name]] <- result_auc
      SVM_result_Youden[[name]] <- result_Youden
      SVM_result_kappa[[name]] <- result_kappa
      SVM_result_model[[name]] <- model
      SVM_result_summary[[name]] <- all_result
    }}
  SVM_result <- list(F1 = SVM_result_F1,
                     ACC = SVM_result_ACC,
                     RECALL = SVM_result_Recall,
                     NPV = SVM_result_NPV,
                     PPV = SVM_result_PPV,
                     Sensitive=SVM_result_Sen,
                     Specific=SVM_result_Spe,
                     AUC=SVM_result_auc,
                     Youden=SVM_result_Youden,
                     Kappa=SVM_result_kappa,
                     MODEL =SVM_result_model,
                     SUMMARY =SVM_result_summary)
  return(SVM_result)

}
