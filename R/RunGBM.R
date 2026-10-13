#' Run Gradient Boosting Machine Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct gradient
#'   boosting machine models.
#' @param fold Number of folds used to estimate the training sample size
#'   for calculating the bag fraction.
#' @param cutoff_method Method used to determine the classification cutoff.
#'   One of \code{"fixed"}, \code{"youden"}, or
#'   \code{"closest.topleft"}.
#' @param cutoff Fixed probability cutoff used when
#'   \code{cutoff_method = "fixed"}.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted gradient
#'   boosting machine models, and prediction results.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- RunGBM(
#'   Object = Object,
#'   gene = gene,
#'   fold = 10,
#'   cutoff_method = "youden",
#'   cutoff = 0.5,
#'   seed = 123
#' )
#' }
RunGBM <- function(Object,
                   gene,
                   fold=10,
                   cutoff_method = c("fixed","youden","closest.topleft"),
                   cutoff=0.5,
                   seed=123){
  #######################################
  GBM_result <- list()
  GBM_result_F1 <- list()
  GBM_result_ACC <- list()
  GBM_result_Recall <- list()
  GBM_result_NPV <- list()
  GBM_result_PPV <- list()
  GBM_result_Sen <-list()
  GBM_result_Spe <- list()
  GBM_result_auc <- list()
  GBM_result_summary <- list()
  GBM_result_Youden <-list()
  GBM_result_kappa <- list()
  GBM_result_model<- list()
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
        "Model: Grandient Boosting Machine","\n",
        "Gene Set:", selectionname, "\n",
        "Gene Number:", length(selection), "\n",
        "====================================\n")
    train_exp <- train_data[,selection,drop=F]
    train_exp$labels <- train_label$Label
    set.seed(seed)
    n_cv_train <- floor(nrow(train_exp) * (fold - 1) / fold)
    min_bag <- 22 / n_cv_train
    bag.fraction <- min(1,max(0.5, min_bag + 0.05))
    model <- gbm::gbm(labels~., data = train_exp,
                      distribution = "bernoulli",
                      bag.fraction = bag.fraction ,
                      n.minobsinnode = 10)
    train_result <- as.data.frame(gbm::predict.gbm(model, train_data[,selection,drop=F],type = "response"))
    train_result$real_label <- factor(ifelse(train_exp$labels==1, "positive", "negative"))
    names(train_result) <- c("predict_p","real_label")
    roc_obj <- pROC::roc(response = train_result$real_label,predictor = train_result$predict_p, levels = c("negative", "positive"),direction = "<")
    Youden_cut <-  pROC::coords(roc_obj,x = "best",best.method = "youden",ret = c("threshold", "sensitivity", "specificity"))
    AUC_cut <-  pROC::coords(roc_obj,x = "best",best.method = "closest.topleft",ret = c("threshold", "sensitivity", "specificity"))
    ###############cutoff:#Youden#or#AUC###########################################################################################################################
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
      tresult <- as.data.frame(predict(model, as.data.frame(expdata),type = "response"))
      tresult$type <- factor(ifelse(tresult[,1] > cutoff_use, "positive", "negative") )
      #zsc
      colnames(tresult) <- c("predict_p", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      rownames(tresult)<-comd
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
    name=paste0(selectionname," + GBM(cutoff:",round(cutoff_use, 2), ")")
    GBM_result_F1[[name]] <- result_F1
    GBM_result_ACC[[name]] <- result_ACC
    GBM_result_Recall[[name]] <- result_Recall
    GBM_result_NPV[[name]] <- result_NPV
    GBM_result_PPV[[name]] <- result_PPV
    GBM_result_Sen[[name]] <- result_Sen
    GBM_result_Spe[[name]] <- result_Spe
    GBM_result_auc[[name]] <- result_auc
    GBM_result_Youden[[name]] <- result_Youden
    GBM_result_kappa[[name]] <- result_kappa
    GBM_result_model[[name]] <- model
    GBM_result_summary[[name]] <- all_result
  }
  GBM_result <- list(F1 = GBM_result_F1,
                     ACC = GBM_result_ACC,
                     RECALL = GBM_result_Recall,
                     NPV = GBM_result_NPV,
                     PPV = GBM_result_PPV,
                     Sensitive=GBM_result_Sen,
                     Specific=GBM_result_Spe,
                     AUC=GBM_result_auc,
                     Youden=GBM_result_Youden,
                     Kappa=GBM_result_kappa,
                     MODEL =GBM_result_model,
                     SUMMARY =GBM_result_summary)
  return(GBM_result)
}
