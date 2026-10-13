#' Run Random Forest SRC Classification Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct random forest
#'   classification models.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted
#'   randomForestSRC models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunRSF(
#'   Object = Object,
#'   gene = gene,
#'   seed = 123
#' )
#' }
RunRSF <- function(Object,
                   gene,
                   seed=123){
  #######################################
  RSF_result <- list()
  RSF_result_F1 <- list()
  RSF_result_ACC <- list()
  RSF_result_Recall <- list()
  RSF_result_NPV <- list()
  RSF_result_PPV <- list()
  RSF_result_Sen <-list()
  RSF_result_Spe <- list()
  RSF_result_auc <- list()
  RSF_result_summary <- list()
  RSF_result_Youden <-list()
  RSF_result_kappa <- list()
  RSF_result_model<- list()
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
          "Model: Fast Unified Random Forests(RSF)","\n",
          "Gene Set:", selectionname, "\n",
          "Gene Number:", length(selection), "\n",
          "mtry=", mnumberi,names(which(mnum==mnumberi)),"p","\n",
          "====================================\n")
      train_exp <- train_data[,selection,drop=F]
      train_exp$labels <- factor(train_label$Label,levels = c(0,1))
      set.seed(seed)
      model <- randomForestSRC::rfsrc(formula = labels~.,
                                      data = train_exp,
                                      ntree = 1000,
                                      mtry=mnumberi,
                                      nodesize = 5,
                                      importance = T,
                                      proximity = T,
                                      forest = T)
      optionTrees <- which.min(model$err.rate[,1])
      set.seed(seed)
      model <- randomForestSRC::rfsrc(labels~.,
                                      data = train_exp,importance = T,
                                      proximity = T,
                                      forest = T,
                                      ntree=optionTrees, mtry=mnumberi)
      train_pred <- predict(model,newdata = train_exp)
      train_result <- data.frame(
        predict_p = train_pred$predicted[, "1"],
        real_label = factor(train_exp$labels,levels = c("0", "1"),labels = c("negative", "positive")))
      train_result$predict_result <- train_pred$class
      train_result$predict_result <- ifelse(train_result$predict_result%in%1,"positive","negative")

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
        tresult_pred <- predict(model, as.data.frame(expdata))
        tresult <- data.frame(
          predict_p = tresult_pred$predicted[, "1"],
          real_label = factor(
            ifelse(real_label_numeric == 1,"positive","negative"),
            levels = c("negative", "positive")),
          predict_result = factor(
            ifelse(tresult_pred$class ==1 ,
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
      name=paste0(selectionname," + RSF(mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")
      RSF_result_F1[[name]] <- result_F1
      RSF_result_ACC[[name]] <- result_ACC
      RSF_result_Recall[[name]] <- result_Recall
      RSF_result_NPV[[name]] <- result_NPV
      RSF_result_PPV[[name]] <- result_PPV
      RSF_result_Sen[[name]] <- result_Sen
      RSF_result_Spe[[name]] <- result_Spe
      RSF_result_auc[[name]] <- result_auc
      RSF_result_Youden[[name]] <- result_Youden
      RSF_result_kappa[[name]] <- result_kappa
      RSF_result_model[[name]] <- model
      RSF_result_summary[[name]] <- all_result
    }}
  RSF_result <- list(
    F1 = RSF_result_F1,
    ACC = RSF_result_ACC,
    RECALL = RSF_result_Recall,
    NPV = RSF_result_NPV,
    PPV = RSF_result_PPV,
    Sensitive=RSF_result_Sen,
    Specific=RSF_result_Spe,
    AUC=RSF_result_auc,
    Youden=RSF_result_Youden,
    Kappa=RSF_result_kappa,
    MODEL =RSF_result_model,
    SUMMARY =RSF_result_summary)
  return(RSF_result)
}
