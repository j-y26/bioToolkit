# Calculate the correlation between gene expression and another variable using a linear model

This function calculates the correlation between gene expression and
another continuous variable, most commonly another gene. This function
is used to identify the list of genes in which the expression is
correlated with another variable/gene in a dose- dependent manner. This
function uses a linear model to calculate the correlation, with the
query gene as the independent variable and all other genes as the
dependent variable. The significance of the correlation is calculated
using the p-value of the linear model and is adjusted for multiple
testing using a user-defined method. The slope of the linear model is
used as the correlation coefficient, with R-squared as the measure of
goodness of fit. The strength of the correlation is calculated using a
user-defined method (e.g., Pearson, Spearman, Kendall).

## Usage

``` r
calc_expr_corr(
  gene_expr,
  query_var,
  metadata = NULL,
  vars_to_regress = NULL,
  scale_expr = TRUE,
  plot_model = FALSE,
  base_dir = getwd(),
  corr_method = "pearson",
  p_adjust_method = "BH",
  file_type = "png",
  dpi = 1000,
  plot_group_var = NULL,
  switch_axes = FALSE,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- gene_expr:

  A data frame containing the gene expression data. The data frame
  should have the genes as rows and samples as columns. The gene
  expression data should be normalized (e.g. CPM) and consist of only
  numeric values. All genes in this data frame, except the query gene,
  will be used as the dependent variable in the linear model. However,
  the query_var can also be specified in the `metadata` data frame.
  Note, only numeric values are allowed in the gene expression data.

- query_var:

  A character string specifying the query gene for which the correlation
  with other genes will be calculated. The query gene should be present
  in the gene expression data frame, and it must be a row name in the
  `gene_expr` data frame. In general, the query gene is the gene of
  interest and should be a continuous variable specified in the
  `gene_expr` data frame. However, in other cases, the query gene can be
  a categorical/ordinal variable specified in the `metadata` data frame.
  In this case, any ordinal variable should be converted into a numeric
  variable before using this function. Categorical variables should be
  converted into a factor and have the base level set as the reference
  level.

- metadata:

  A data frame containing the metadata information for the samples. The
  data frame should have samples as rows and metadata variables as
  columns. The metadata information can contain the `query_var` and/or
  `vars_to_regress` variables that are used as covariates in the linear
  model. The metadata information can be used to adjust the gene
  expression data for confounding variables. The metadata can contain
  both categorical and continuous variables. ordinal variables should be
  converted into numeric variables before using this function.
  Categorical variables should be converted into a factor and have the
  base level set as the reference level. Default value is `NULL`. Note,
  row names of the `metadata` data frame should match the column names
  of the `gene_expr` data frame.

- vars_to_regress:

  A character vector specifying the variables to regress out from the
  linear model. The default value is `NULL`, which means no variables
  will be regressed out. The variables can present in either the
  `gene_expr` or `metadata` data frame. However, if the variable are not
  continuous, they should only exist in the `metadata` data frame.

- scale_expr:

  A logical value indicating whether to scale the gene expression data
  before calculating the correlation. If `TRUE`, the gene expression
  data will be scaled using the `scale` function to a mean of 0 and a
  standard deviation of 1 (z-score). Expression of each gene will be
  scaled independently across samples. The default value is `TRUE`. If
  false, the original gene expression data will be used for
  calculations.

- plot_model:

  A logical value indicating whether to plot the linear model for each
  gene. The default value is `FALSE`. If true, a plot will be generated
  for each gene showing the linear model with the query gene on the
  x-axis and the gene expression on the y-axis. The plot will also show
  the linear model coefficients, R-squared value, correlation
  coefficient, and FDR-adjusted p-value. The plot will be a ggplot
  object and can be further customized using the ggplot2 functions. The
  plot will be saved in the `model_plots_positive` and
  `model_plots_negative` directories in the `base_dir` directory.

- base_dir:

  A character string specifying the base directory to save the model
  results and plots. The default value is the current working directory.
  The model results will be saved as a CSV file named
  `query_var_correlation_results.csv`.

- corr_method:

  A character string specifying the method to calculate the correlation.
  The available methods include "pearson", "spearman", "kendall". The
  default method is "pearson".

- p_adjust_method:

  A character string specifying the method to adjust the p-values for
  multiple testing. The available methods include "BH"
  (Benjamini-Hochberg), "bonferroni", "holm", "hochberg", "hommel",
  "BY", "fdr", "none". The default method is "BH". See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for more details.

- file_type:

  A character string specifying the file type for the model plots. The
  available file types include "png", "jpeg", "tiff", "pdf", "svg". The
  default file type is "png".

- dpi:

  An integer specifying the resolution of raster plots (e.g. png) in
  dots per inch (DPI). The default value is 1000. The DPI value should
  be greater than 0.

- plot_group_var:

  A character string specifying the variable in the metadata data frame
  to group the samples in the plot.

- switch_axes:

  A logical value indicating whether to switch the axes of the plot. The
  default value is `FALSE`. If `TRUE`, the query gene will be plotted on
  the y-axis and the gene expression on the x-axis.

- alpha:

  The significance level for the adjusted p-values. The default value is
  0.05.

- verbose:

  A logical value indicating whether to print the progress messages. The
  default value is TRUE.

## Value

A list of four elements:

- model_list:

  A list of linear models for each gene, with the query gene as the
  independent variable and the gene expression as the dependent
  variable.

- model_results:

  A data frame containing the correlation results, including the
  correlation coefficient, p-value, adjusted p-value, and R-squared
  value for each gene.

- p_volcano:

  A ggplot object of the volcano plot showing the correlation results.

- pos_plots:

  A list containing the ggplot objects for each gene showing the linear
  model plot. The plots will be saved automatically. Only available if
  `plot_model` is `TRUE`.

- neg_plots:

  A list containing the ggplot objects for each gene showing the linear
  model plot. The plots will be saved automatically. Only available if
  `plot_model` is `TRUE`.
