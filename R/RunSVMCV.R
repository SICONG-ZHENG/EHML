#' Run Support Vector Machine Models with Cross-Validation
#'
#' @description
#' Trains and evaluates support vector machine (SVM) models
#' using cross-validation and multiple kernel functions.
#'
#' The function supports linear, radial basis function,
#' and polynomial kernels. Models are constructed separately
#' for each supplied feature set.
#'
#' @param Object A list containing training data, training labels,
#'   validation datasets, and validation labels.
#'
#' @param gene A named list of character vectors. Each element
#'   specifies the predictor variables used to construct an SVM model.
#'
#' @param kernel_all A character vector specifying the SVM methods
#'   evaluated by `caret`. Available options include `"svmLinear"`,
#'   `"svmRadial"`, and `"svmPoly"`.
#'
#' @param fold A positive integer specifying the number of
#'   cross-validation folds. Defaults to 10.
#'
#' @param seed A positive integer specifying the random seed
#'   for reproducible model training. Defaults to 123.
#'
#' @return A list containing model performance metrics,
#'   fitted SVM models, and sample-level prediction results.
#'
#' @details
#' Support vector machine models are trained using the
#' `caret` framework with cross-validation.
#'
#' Different kernel functions are evaluated to capture
#' linear and nonlinear relationships between predictors
#' and binary outcomes.
#'
#' @export
RunSVMCV<- function(Object,
                    gene,
                    kernel_all =c("svmLinear", "svmRadial", "svmPoly"),
                    fold=10,
                    seed=123){
  #######################################
  SVMCV_result <- list()
  SVMCV_result_F1 <- list()
  SVMCV_result_ACC <- list()
  SVMCV_result_Recall <- list()
  SVMCV_result_NPV <- list()
  SVMCV_result_PPV <- list()
  SVMCV_result_Sen <-list()
  SVMCV_result_Spe <- list()
  SVMCV_result_auc <- list()
  SVMCV_result_summary <- list()
  SVMCV_result_Youden <-list()
  SVMCV_result_kappa <- list()
  SVMCV_result_model<- list()
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
      train_control <- trainControl(method="cv", number=fold)
      model <- caret::train(labels~., data = train_exp,
                            method = kerneli,
                            trControl = train_control)
      model <- model$finalModel
      pred <- kernlab::predict(model, as.data.frame(train_data[,selection,drop=F]))
      prob <-  kernlab::predict(model,as.data.frame(train_data[,selection,drop=F]),type = "decision")
      train_result <- data.frame(real_label=factor(ifelse(train_exp$labels==1, "positive", "negative")),
                                 predict_result=factor(ifelse(pred== 1, "positive", "negative") ),
                                 predict_p = prob[,1])
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
        tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
        tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
        colnames(tresult) <- c("real_label", "predict_result")
        tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
        #zsc
        probt <-  kernlab::predict(model,as.data.frame(expdata),type = "decision")
        tresult$predict_p <-probt[,1]
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
      name=paste0(selectionname," + SVM-CV(kernel:",kerneli,", ",fold,"fold)")
      SVMCV_result_F1[[name]] <- result_F1
      SVMCV_result_ACC[[name]] <- result_ACC
      SVMCV_result_Recall[[name]] <- result_Recall
      SVMCV_result_NPV[[name]] <- result_NPV
      SVMCV_result_PPV[[name]] <- result_PPV
      SVMCV_result_Sen[[name]] <- result_Sen
      SVMCV_result_Spe[[name]] <- result_Spe
      SVMCV_result_auc[[name]] <- result_auc
      SVMCV_result_Youden[[name]] <- result_Youden
      SVMCV_result_kappa[[name]] <- result_kappa
      SVMCV_result_model[[name]] <- model
      SVMCV_result_summary[[name]] <- all_result
    }}
  SVMCV_result <- list(F1 = SVMCV_result_F1,
                       ACC = SVMCV_result_ACC,
                       RECALL = SVMCV_result_Recall,
                       NPV = SVMCV_result_NPV,
                       PPV = SVMCV_result_PPV,
                       Sensitive=SVMCV_result_Sen,
                       Specific=SVMCV_result_Spe,
                       AUC=SVMCV_result_auc,
                       Youden=SVMCV_result_Youden,
                       Kappa=SVMCV_result_kappa,
                       MODEL =SVMCV_result_model,
                       SUMMARY =SVMCV_result_summary)
  return(SVMCV_result)
}
