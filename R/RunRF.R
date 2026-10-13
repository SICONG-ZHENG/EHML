#' Run Random Forest Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct random forest models.
#' @param fold Number of cross-validation folds. Currently not used in the function.
#' @param seed Random seed used to ensure reproducible random forest models.
#'
#' @return A list containing model performance metrics, fitted random forest
#'   models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunRF(
#'   Object = Object,
#'   gene = gene,
#'   seed = 123
#' )
#' }
RunRF<- function(Object,
                 gene,
                 fold=10,
                 seed=123){
  #######################################
  RF_result <- list()
  RF_result_F1 <- list()
  RF_result_ACC <- list()
  RF_result_Recall <- list()
  RF_result_NPV <- list()
  RF_result_PPV <- list()
  RF_result_Sen <-list()
  RF_result_Spe <- list()
  RF_result_auc <- list()
  RF_result_summary <- list()
  RF_result_Youden <-list()
  RF_result_kappa <- list()
  RF_result_model<- list()
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
    mnum <- c(round(quantile(1:length(selection), probs = 0.25)),
              round(quantile(1:length(selection), probs = 0.5)),
              round(quantile(1:length(selection), probs = 0.75)))
    mnum <- mnum[!duplicated(mnum)]
    for (mnumberi in mnum) {
      cat("\n====================================\n",
          "Model: Random forest(RF)","\n",
          "Gene Set:", selectionname, "\n",
          "Gene Number:", length(selection), "\n",
          "mtry=", mnumberi,names(which(mnum==mnumberi)),"p","\n",
          "====================================\n")
      train_exp <- train_data[,selection,drop=F]
      train_exp$labels <- train_label$Label
      train_exp$labels <- factor(train_exp$labels,levels = c(0,1))
      set.seed(seed)
      model <- randomForest::randomForest(labels~., data = train_exp, ntree=1000, mtry=mnumberi)
      optionTrees <- which.min(model$err.rate[,1])  #
      set.seed(123)
      model <- randomForest::randomForest(labels~.,
                                          data = train_exp,
                                          ntree=optionTrees, mtry=mnumberi)
      train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
      train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
      colnames(train_result) <- c("real_label", "predict_result")
      train_result$real_label <- factor(ifelse(train_exp$labels==1, "positive", "negative"))
      #zsc
      train_result$predict_p <-predict(model, as.data.frame(train_exp), type = "prob")[,2]
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
        tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
        tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
        colnames(tresult) <- c("real_label", "predict_result")
        tresult$real_label <- factor(ifelse(real_label_numeric==1, "positive", "negative"))
        #zsc
        tresult$predict_p <-predict(model, as.data.frame(expdata), type = "prob")[,2]
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
      name=paste0(selectionname," + RF(mtry:", mnumberi,", ", names(which(mnum==mnumberi)), "p)")
      RF_result_F1[[name]] <- result_F1
      RF_result_ACC[[name]] <- result_ACC
      RF_result_Recall[[name]] <- result_Recall
      RF_result_NPV[[name]] <- result_NPV
      RF_result_PPV[[name]] <- result_PPV
      RF_result_Sen[[name]] <- result_Sen
      RF_result_Spe[[name]] <- result_Spe
      RF_result_auc[[name]] <- result_auc
      RF_result_Youden[[name]] <- result_Youden
      RF_result_kappa[[name]] <- result_kappa
      RF_result_model[[name]] <- model
      RF_result_summary[[name]] <- all_result
    }}
  RF_result <- list(F1 = RF_result_F1,
                    ACC = RF_result_ACC,
                    RECALL = RF_result_Recall,
                    NPV = RF_result_NPV,
                    PPV = RF_result_PPV,
                    Sensitive=RF_result_Sen,
                    Specific=RF_result_Spe,
                    AUC=RF_result_auc,
                    Youden=RF_result_Youden,
                    Kappa=RF_result_kappa,
                    MODEL =RF_result_model,
                    SUMMARY =RF_result_summary)
  return(RF_result)

}
