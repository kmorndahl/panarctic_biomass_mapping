bias_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  
  bias_impl <- function(truth, estimate) {
    
    mean(estimate-truth)

  }
  
  metric_vec_template(
    metric_impl = bias_impl,
    truth = truth, 
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric",
    ...
  )
  
}

bias <- function(data, ...) {
  UseMethod("bias")
}

bias <- new_numeric_metric(bias, direction = "minimize")

bias.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  numeric_metric_summarizer(
    name = "bias",
    fn = bias_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate), 
    na_rm = na_rm,
    ...
  )
  
}

bias_vec(
  truth = solubility_test$solubility, 
  estimate = solubility_test$prediction
)

bias(solubility_test, truth = solubility, estimate = prediction)
