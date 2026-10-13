#' Run Elastic Net Regression Models with Cross-Validation
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct elastic net
#'   regression models.
#' @param fold Number of folds used for cross-validation.
#' @param alpha Mixing parameter of the elastic net penalty.
#'   Values range from 0 (ridge regression) to 1 (lasso regression).
#' @param cutoff_method Method used to determine the classification cutoff.
#'   One of \code{"fixed"}, \code{"youden"}, or
#'   \code{"closest.topleft"}.
#' @param cutoff Fixed probability cutoff used when
#'   \code{cutoff_method = "fixed"}.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted elastic net
#'   regression models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunEnetCV(
#'   Object = Object,
#'   gene = gene,
#'   fold = 10,
#'   alpha = seq(0.1, 0.9, 0.1),
#'   cutoff_method = "youden",
#'   cutoff = 0.5,
#'   seed = 123
#' )
#' }
RunEnetCV <- function(Object,
                      gene,
                      fold=10,
                      alpha=seq(0.1,0.9,0.1),
                      cutoff_method = c("fixed","youden","closest.topleft"),
                      cutoff=0.5,
                      seed=123){
  #######################################
  EnetCV_result <- list()
  EnetCV_result_F1 <- list()
  EnetCV_result_ACC <- list()
  EnetCV_result_Recall <- list()
  EnetCV_result_NPV <- list()
  EnetCV_result_PPV <- list()
  EnetCV_result_Sen <-list()
  EnetCV_result_Spe <- list()
  EnetCV_result_auc <- list()
  EnetCV_result_summary <- list()
  EnetCV_result_Youden <-list()
  EnetCV_result_kappa <- list()
  EnetCV_result_model<- list()
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
    for (alphai in alpha) {
      selection <- gene[[i]]
      selectionname <- names(gene)[i]
      cat("\n====================================\n",
          "Model: Elastic Net regression-CV","\n",
          "Gene Set:", selectionname, "\n",
          "Gene Number:", length(selection), "\n",
          "Alpha:",alphai, "\n",
          "====================================\n")
      train_exp <- as.matrix(train_data[,selection,drop=F])
      labels <- factor(train_label$Label)
      set.seed(seed)
      cvfit <- glmnet::cv.glmnet(train_exp, labels,family = "binomial", nlambda=100,
                                 alpha=alphai,nfolds = fold)
      train_result <- as.data.frame(predict(cvfit,
                                            train_exp,
                                            s=cvfit$lambda.min, type="response"))
      train_result$real_label <- factor(ifelse(labels==1, "positive", "negative"))
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
        tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
        predict_p <-predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response")
        tresult <- data.frame(
          predict_p = predict_p[,1],
          real_label = factor(
            ifelse(real_label_numeric == 1,"positive","negative"),
            levels = c("negative", "positive")),
          predict_result = factor(
            ifelse(predict_p[,1] > cutoff_use,
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
      name=paste0(selectionname," + ENR-CV(alpha:",alphai,", cutoff:",  round(cutoff_use,2),", ", fold, "fold)")
      EnetCV_result_F1[[name]] <- result_F1
      EnetCV_result_ACC[[name]] <- result_ACC
      EnetCV_result_Recall[[name]] <- result_Recall
      EnetCV_result_NPV[[name]] <- result_NPV
      EnetCV_result_PPV[[name]] <- result_PPV
      EnetCV_result_Sen[[name]] <- result_Sen
      EnetCV_result_Spe[[name]] <- result_Spe
      EnetCV_result_auc[[name]] <- result_auc
      EnetCV_result_Youden[[name]] <- result_Youden
      EnetCV_result_kappa[[name]] <- result_kappa
      EnetCV_result_model[[name]] <- cvfit
      EnetCV_result_summary[[name]] <- all_result
    }}
  EnetCV_result <- list(F1 = EnetCV_result_F1,
                        ACC = EnetCV_result_ACC,
                        RECALL = EnetCV_result_Recall,
                        NPV = EnetCV_result_NPV,
                        PPV = EnetCV_result_PPV,
                        Sensitive=EnetCV_result_Sen,
                        Specific=EnetCV_result_Spe,
                        AUC=EnetCV_result_auc,
                        Youden=EnetCV_result_Youden,
                        Kappa=EnetCV_result_kappa,
                        MODEL =EnetCV_result_model,
                        SUMMARY =EnetCV_result_summary)
  return(EnetCV_result)
}
