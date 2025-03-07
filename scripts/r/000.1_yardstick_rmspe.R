library(rlang)

rmspe_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  
  rmspe_impl <- function(truth, estimate) {
    
    r = ((truth - estimate)/truth)^2
    sqrt(mean(r[is.finite(r)]))*100
    
  }
  
  metric_vec_template(
    metric_impl = rmspe_impl,
    truth = truth, 
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric",
    ...
  )
  
}

rmspe <- function(data, ...) {
  UseMethod("rmspe")
}

rmspe <- new_numeric_metric(rmspe, direction = "minimize")

rmspe.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  numeric_metric_summarizer(
    name = "rmspe",
    fn = rmspe_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate), 
    na_rm = na_rm,
    ...
  )
  
}

rmspe_vec(
  truth = solubility_test$solubility, 
  estimate = solubility_test$prediction
)

rmspe(solubility_test, truth = solubility, estimate = prediction)
