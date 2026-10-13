

#' Split Data into Training and Validation Datasets
#' @description
#' Splits an expression or feature matrix into one training dataset and
#' multiple validation datasets according to dataset information provided
#' in a sample annotation table. Binary outcome labels are converted to
#' numeric values, where the class specified by \code{Positive} is coded
#' as 1 and the other class is coded as 0.
#'
#' @param Data A matrix or data frame containing samples in rows and
#'   variables in columns. Row names must be sample identifiers.
#' @param Group A data frame containing sample annotation information.
#'   Row names must correspond to the sample identifiers in \code{Data}.
#' @param Data_name A character string specifying the column in
#'   \code{Group} that contains dataset names, such as \code{"Data"}.
#' @param Valiabel_name A character string specifying the column in
#'   \code{Group} that contains the binary outcome label.
#' @param Train_names A character vector specifying which dataset or
#'   datasets should be used as the training set. Multiple dataset names
#'   can be supplied and will be combined into one training dataset.
#' @param Positive The value representing the positive class. The specified
#'   class is converted to 1, and the other class is converted to 0.
#'   Default is 1.
#'
#' @returns A list containing:
#'   \describe{
#'   \item{\code{train_data}}{A data frame containing the training samples.}
#'   \item{\code{train_label}}{A one-column data frame containing the
#'     binary training labels.}
#'   \item{\code{test_list}}{A named list of validation datasets.}
#'   \item{\code{test_label}}{A named list containing the corresponding
#'     binary labels for each validation dataset.}
#' }
#'
#' @details
#' Both \code{Data} and \code{Group} must contain unique sample identifiers
#' as row names. Only binary classification outcomes are supported.
#'
#' Labels may be numeric, integer, factor, character, or logical. Factor
#' labels are first converted to character values before the positive class
#' is identified.
#'
#' @export
#'
#' @examples
#'
#' data(iris)
#'
#' rownames(iris) <- paste0("Iris", seq_len(nrow(iris)))
#'
#' Data <- iris[iris$Species != "setosa", 1:4]
#' Group <- iris[iris$Species != "setosa", "Species", drop = FALSE]
#'
#' Group$Data <- rep(
#'   c("Train", "Test1", "Test2", "Test3"),
#'   length.out = nrow(Group)
#' )
#'
#' result <- DataChange(
#'   Data = Data,
#'   Group = Group,
#'   Data_name = "Data",
#'   Valiabel_name = "Species",
#'   Train_names = "Train",
#'   Positive = "versicolor"
#' )
#'
#' names(result)
#' names(result$test_list)
#' table(result$train_label$Label)
#'
DataChange <- function(Data,
                       Group,
                       Data_name,
                       Valiabel_name,
                       Train_names,
                       Positive = 1){
  dataset_vector <- as.character(Group[, Data_name])
  label_vector <- Group[, Valiabel_name]
  if(is.factor(label_vector)){
    label_vector <- as.character(label_vector)
  }
  if(is.character(label_vector)){

    label_vector <- ifelse(label_vector == Positive,1,0)
  } else {
    u <- sort(unique(label_vector))
    if(length(u) != 2){
      stop("Only binary classification is supported.")
    }
    if(all(u %in% c(0,1))){
      label_vector <- as.numeric(label_vector)
    }else{
      label_vector <- ifelse(
        label_vector == Positive, 1,0)}}

  ################################################
  # Train
  ################################################
  train_id <- rownames(Group)[dataset_vector %in% Train_names]
  train_data <- Data[train_id,,drop = FALSE]
  train_label <- data.frame(Label = label_vector[dataset_vector %in% Train_names])
  rownames(train_label) <- train_id

  ################################################
  # Test
  ################################################
  test_names <- setdiff(
    unique(dataset_vector),
    Train_names
  )
  test_list <- list()
  test_label <- list()
  for(i in seq_along(test_names)){
    db <- test_names[i]
    sample_id <- rownames(Group)[dataset_vector == db]
    test_list[[db]] <- Data[sample_id,,drop = FALSE]
    tmp_lab <- data.frame(Label = label_vector[dataset_vector == db])
    rownames(tmp_lab) <- sample_id
    test_label[[db]] <- tmp_lab}
  ################################################
  # return
  ################################################
  return(list(train_data = train_data,train_label = train_label,test_list = test_list,test_label = test_label
  )
  )
}
