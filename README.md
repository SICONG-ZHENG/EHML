# EHML: Easy Hyperparameter Optimization for Machine Learning

**EHML** is an R-based framework designed to facilitate machine learning model development, hyperparameter optimization, feature selection, model evaluation, and performance visualization.

The framework integrates multiple machine learning algorithms into a unified workflow for binary classification tasks, with a particular focus on biomedical and clinical research.

## Key Features

- **Data preparation:** Organization of training and validation datasets for binary classification.
- **Feature selection:** Integration of Elastic Net, Lasso, stepwise logistic regression, GLMBoost, sparse PLS-GLM, random forest, GBM, and Boruta-based methods.
- **Machine learning:** Implementation of statistical and machine learning algorithms, including logistic regression, discriminant analysis, support vector machines, random forests, gradient boosting, XGBoost, and neural networks.
- **Hyperparameter optimization:** Evaluation of multiple model configurations and parameter combinations.
- **Model evaluation:** Assessment using AUC, accuracy, F1 score, recall, sensitivity, specificity, positive predictive value, negative predictive value, Youden index, and Cohen's kappa.
- **Visualization:** Visualization and comparison of predictive performance across training and validation datasets, including an Archimedean spiral heatmap.

## Installation

EHML is currently under development and is hosted on GitHub.

```r
install.packages("remotes")

remotes::install_github("SICONG-ZHENG/EHML")
```

After installation:

```r
library(EHML)
```

**Note:** The repository is currently private. Installation requires authorized GitHub access.

## Requirements

EHML is developed in R and uses established R packages for its machine learning algorithms.

Some modules have additional requirements:

- **PLS-GLM:** `plsRglm` version 1.5.1 is recommended.
- **XGBoost:** `xgboost` version 1.5.2.1 is recommended.
- **Neural networks:** A separately configured Python environment with compatible versions of Keras and TensorFlow is required.

Python should be configured before using the neural network module.

## General Workflow

The EHML analysis workflow consists of the following stages:

1. Data preparation and organization of training and validation datasets.
2. Feature selection using multiple statistical and machine learning approaches.
3. Hyperparameter optimization and model construction.
4. Performance evaluation across training and validation datasets.
5. Model comparison and performance visualization.

## Documentation

Function-level documentation is available through the R help system.

For example:

```r
features
RunLR
RunRF
RunNN
MLplot
```

## Development Status

EHML is currently under active development. The software and associated methodology are being prepared for scientific publication.

Compatibility and reproducibility testing across different R and Python environments is ongoing.

## License

EHML is distributed under the MIT License.

## Author

Sicong Zheng

## Citation

A manuscript describing the EHML framework is currently in preparation. Citation information will be updated following publication.
