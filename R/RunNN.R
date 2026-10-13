#' Run Neural Network (MLP) Models
#'
#' @param Object A list containing training data, training labels,
#'   test datasets, and test labels.
#' @param gene A named list of gene sets used to construct neural
#'   network (MLP) models.
#' @param epochselecti Number of training epochs.
#' @param learningrate A numeric vector of learning rates to evaluate.
#' @param dropoutrate A numeric vector of dropout rates to evaluate.
#' @param envname Name of the Python virtual environment containing
#'   TensorFlow and Keras.
#' @param cutoff_method Method used to determine the classification cutoff.
#'   One of \code{"fixed"}, \code{"youden"}, or
#'   \code{"closest.topleft"}.
#' @param cutoff Fixed probability cutoff used when
#'   \code{cutoff_method = "fixed"}.
#' @param seed Random seed used to ensure reproducible model training.
#'
#' @return A list containing model performance metrics, fitted neural
#'   network (MLP) models, and prediction results.
#' @export
#' @importFrom keras keras_model_sequential
#' @importFrom keras layer_dense layer_dropout
#' @importFrom keras compile optimizer_adam fit
#' @importFrom dplyr %>%
#' @importFrom stats predict
#' @details
#' Saved models are exported in the native Keras (.keras) format.
#'
#' @examples
#' \dontrun{
#' result <- RunNN(
#'   Object = Object,
#'   gene = gene,
#'   epochselecti = 50,
#'   learningrate = c(0.001, 0.005, 0.01, 0.05),
#'   dropoutrate = c(0.25, 0.5, 0.75),
#'   envname = "r-tensorflow",
#'   cutoff_method = "youden",
#'   cutoff = 0.5,
#'   seed = 123
#' )
#' }
RunNN <- function(Object,
                  gene,
                  epochselecti = 50,
                  learningrate = c(0.001, 0.005, 0.01, 0.05),
                  dropoutrate = c(0.25, 0.5, 0.75),
                  envname="r-tensorflow",
                  cutoff_method = c("fixed","youden","closest.topleft"),
                  cutoff=0.5,
                  seed=123){
  #######################################
  NN_result <- list()
  NN_result_F1 <- list()
  NN_result_ACC <- list()
  NN_result_Recall <- list()
  NN_result_NPV <- list()
  NN_result_PPV <- list()
  NN_result_Sen <-list()
  NN_result_Spe <- list()
  NN_result_auc <- list()
  NN_result_summary <- list()
  NN_result_Youden <-list()
  NN_result_kappa <- list()
  NN_result_model<- list()
  #########
  if (!dir.exists("Model for Neural Network (MLP)")) {
    dir.create("Model for Neural Network (MLP)")
  }
  ########
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
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("Please install the 'keras' package.")
  }

  if (!requireNamespace("tensorflow", quietly = TRUE)) {
    stop("Please install the 'tensorflow' package.")
  }

  reticulate::use_virtualenv(envname, required = TRUE)
  #####
  a <- 2:floor(log2(nrow(train_data))-1)
  batch_size <- 2^a
  #########################################
  Sys.setenv(TF_DETERMINISTIC_OPS = "1")
  for (i in 1:length(gene)) {
    for (batch_sizei in batch_size) {
      for (learningratei in learningrate) {
        for (dropoutratei in dropoutrate) {

          set.seed(seed)
          tensorflow::tf$keras$utils$set_random_seed(as.integer(seed))
          try(tensorflow::tf$config$experimental$enable_op_determinism(),silent = TRUE)
          mlpmodel <- keras_model_sequential()

          selection <- gene[[i]]
          selectionname <- names(gene)[i]
          cat("\n====================================\n",
              "Model: Neural Network (MLP)","\n",
              "Gene Set:", names(gene)[i], "\n",
              "Gene Number:", length(selection), "\n",
              "Batch Size:", batch_sizei, "\n",
              "Learning Rate:", learningratei, "\n",
              "Dropout Rate:", dropoutratei, "\n",
              "====================================\n")

          train_exp <- train_data[,selection,drop=F]
          mlpmodel %>%
            layer_dense(units = 32, activation = "relu", input_shape = ncol(train_exp)) %>%
            layer_dropout(rate = dropoutratei)  #
          mlpmodel %>%
            layer_dense(units = 16, activation = "relu") %>%
            layer_dropout(rate = dropoutratei)
          mlpmodel %>%
            layer_dense(units = 8, activation = "relu") %>%
            layer_dropout(rate = dropoutratei)
          #
          mlpmodel %>%
            layer_dense(units = 1, activation = "sigmoid")
          #
          mlpmodel %>% compile(
            loss = "binary_crossentropy",
            optimizer = optimizer_adam(learning_rate = learningratei),
            metrics = c("accuracy"))
          #
          train_exp <- as.matrix(train_exp)
          history <- mlpmodel %>% fit(
            x = train_exp,
            y = train_label$Label,
            epochs = epochselecti,
            batch_size = batch_sizei)

          summary(mlpmodel)

          train_result <-  as.data.frame(mlpmodel %>% predict(train_exp))
          colnames(train_result) <- "predict_p"
          train_result$real_label <- factor(ifelse(train_label$Label==1, "positive", "negative"),levels = c("negative", "positive"))
          roc_obj <- pROC::roc(response = train_result$real_label,predictor = train_result$predict_p, levels = c("negative", "positive"),direction = "<")
          Youden_cut <-  pROC::coords(roc_obj,x = "best",best.method = "youden",ret = c("threshold", "sensitivity", "specificity"))
          AUC_cut <-  pROC::coords(roc_obj,x = "best",best.method = "closest.topleft",ret = c("threshold", "sensitivity", "specificity"))
          #########################################################################################################################################
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
          train_result$type <- factor(ifelse(train_result[,1] >= cutoff_use, "positive", "negative") )
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
            expdata <- as.matrix(expdata)
            real_label_numeric <- labelsdata[comd,"Label"]
            predict_p <-  as.numeric( mlpmodel %>% predict(expdata))
            tresult <- data.frame(
              predict_p = predict_p,
              real_label = factor(
                ifelse(real_label_numeric == 1,"positive","negative"),
                levels = c("negative", "positive")),
              predict_result = factor(
                ifelse(predict_p > cutoff_use,
                       "positive",
                       "negative"),
                levels = c("negative", "positive")),
              row.names = comd)
            return(tresult)
          },data = test_explist ,labelsdata=test_label)

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
          name=paste0(gsub("[^A-Za-z0-9]", "_", selectionname),"_NN-MLP_", "lr_", learningratei, "_bs_", batch_sizei, "_ep_", epochselecti, "_dropout_", dropoutratei,"_cutoff_",round(cutoff_use, 2))
          NN_result_F1[[name]] <- result_F1
          NN_result_ACC[[name]] <- result_ACC
          NN_result_Recall[[name]] <- result_Recall
          NN_result_NPV[[name]] <- result_NPV
          NN_result_PPV[[name]] <- result_PPV
          NN_result_Sen[[name]] <- result_Sen
          NN_result_Spe[[name]] <- result_Spe
          NN_result_auc[[name]] <- result_auc
          NN_result_Youden[[name]] <- result_Youden
          NN_result_kappa[[name]] <- result_kappa
          NN_result_model[[name]] <- mlpmodel
          NN_result_summary[[name]] <- all_result
          model_path <- paste0("Model for Neural Network (MLP)/",name,".keras")
          mlpmodel$save(model_path)
        }
      }
    }
  }
  NN_result <- list(
    F1 = NN_result_F1,
    ACC = NN_result_ACC,
    RECALL = NN_result_Recall,
    NPV = NN_result_NPV,
    PPV = NN_result_PPV,
    Sensitive=NN_result_Sen,
    Specific=NN_result_Spe,
    AUC=NN_result_auc,
    Youden=NN_result_Youden,
    Kappa=NN_result_kappa,
    MODEL =NN_result_model,
    SUMMARY =NN_result_summary)
  return(NN_result)
}
