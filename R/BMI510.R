
# bmi510.R
#' rando
#' tests whether x is an atomic vector or dataframe-like object
#' @param x
#' @param n, replace
#'
#' @return n samples or n rows
#' @export
#'
#' @examples
rando = function(x, n=1, replace=T){
  if (is.atomic(x)) {
    return(sample(x, n, replace=replace))
  } else if (is.data.frame(x)) {
    return(x[sample(nrow(x), n, replace=replace), ])
  }
}


#' is_min
#' where x equals its maximum value
#' @param x
#' @param na.rm
#'
#' @return a logical with TRUE
#' @export
#'
#' @examples
is_min = function(x, na.rm = T) {
  if (na.rm) {
    x = na.omit(x)
  }
  min_x = min(x)
  return(x == min_x)
}

#' is_max
#' where x equals its maximum value
#' @param x
#' @param na.rm
#'
#' @return a logical with TRUE
#' @export
#'
#' @examples

is_max = function(x, na.rm = T) {
  if (na.rm) {
    x = na.omit(x)
  }
  max_x = max(x)
  return(x == max_x)
}

#' rep_mat
#' make the original matrix become a matrix created by replicating the rows (or columns) M (N) times
#' @param x
#' @param M,N
#'
#' @return a matrix created by replicating the rows (or columns) M (N) times
#' @export
#'
#' @examples

rep_mat = function(x, M=1, N=1){
  if (N > 1) {
    x = t(x)
  }
  rep_x = matrix(NA, nrow=nrow(x)*M, ncol=ncol(x)*N)
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      rep_x[(i-1)*M+1:i*M, (j-1)*N+1:j*N] = x[i,j]
    }
  }
  if (N > 1) {
    rep_x = t(rep_x)
  }
  return(rep_x)
}


#' classes
#' tell people the class of each variable in the vector
#' @param x
#' @param
#'
#' @return a character vector containing the classes of each variable in a tibble x
#' @export
#'
#' @examples

classes = function(x) {
  if (!inherits(x, "tbl_df") && !inherits(x, "data.frame")) {
    stop("x must be a tibble or data frame")
  }
  return(sapply(x, class))
}

#' df_scale
#' scale the tibble x
#' @param x
#' @param center,scale
#'
#' @return a tibble x in which the numeric variables have been scaled with scale
#' @export
#'
#' @examples

df_scale = function(x, center = T, scale = T) {
  num_cols = sapply(x, is.numeric)
  if (any(num_cols)) {
    x[num_cols] = scale(x[num_cols], center = center, scale = scale)
  }
  attr(x, "class") = classes(x)
  return(x)
}


#' log_likelihood_norm
#' calculate the log-likelihood of a sample x under the normal distribution
#' @param x
#' @param mean,sd
#'
#' @return the log-likelihood of a sample x under the normal distribution
#' @export
#'
#' @examples

log_likelihood_norm = function(x, mean, sd) {
  -sum(dnorm(x, mean, sd, log = TRUE))
}

#' log_likelihood_unif
#' calculate the log-likelihood of a sample x under the uniform distribution
#' @param x
#' @param min,max
#'
#' @return the log-likelihood of a sample x under the uniform distribution
#' @export
#'
#' @examples

log_likelihood_unif = function(x, min, max) {
  -sum(dunif(x, min, max, log = TRUE))
}

#' log_likelihood_chisq
#' calculate the log-likelihood of a sample x under the chi-squared
#' @param x
#' @param df
#'
#' @return the log-likelihood of a sample x under the chi-squared
#' @export
#'
#' @examples

log_likelihood_chisq = function(x, df) {
  -sum(dchisq(x, df, log = TRUE))
}

#' log_likelihood_f
#' calculate the log-likelihood of a sample x under the f
#' @param x
#' @param df1,df2
#'
#' @return the log-likelihood of a sample x under the f
#' @export
#'
#' @examples

log_likelihood_f = function(x, df1, df2) {
  -sum(df1/2*log(df1/2) + df2/2*log(df2/2) + (df1/2-1)*log(x) - ((df1+df2)/2)*log(df2+df1*x/df2) - lgamma(df1/2) - lgamma(df2/2) + lgamma((df1+df2)/2))
}

#' log_likelihood_t
#' calculate the log-likelihood of a sample x under the t densities
#' @param x
#' @param df
#'
#' @return the log-likelihood of a sample x under the t densities
#' @export
#'
#' @examples

log_likelihood_t = function(x, df) {
  -sum(dt(x, df, log = TRUE))
}

#' sensitivity
#' Calculate sensitivity based on comparing predicted and training labels
#' @param pred
#' @param truth
#'
#' @return sensitivity
#' @export
#'
#' @examples

sensitivity = function(pred, truth) {
  sum(pred == 1 & truth == 1) / sum(truth == 1)
}

#' specificity
#' Calculate specificity based on comparing predicted and training labels
#' @param pred
#' @param truth
#'
#' @return specificity
#' @export
#'
#' @examples

specificity = function(pred, truth) {
  sum(pred == 0 & truth == 0) / sum(truth == 0)
}

#' precision
#' Calculate precision based on comparing predicted and training labels
#' @param pred
#' @param truth
#'
#' @return precision
#' @export
#'
#' @examples

precision = function(pred, truth) {
  sum(pred == 1 & truth == 1) / sum(pred == 1)
}

#' recall
#' Calculate recall based on comparing predicted and training labels
#' @param pred
#' @param truth
#'
#' @return recall
#' @export
#'
#' @examples

recall = function(pred, truth) {
  sum(pred == 1 & truth == 1) / sum(truth == 1)
}

#' accuracy
#' Calculate accuracy based on comparing predicted and training labels
#' @param pred
#' @param truth
#'
#' @return accuracy
#' @export
#'
#' @examples

accuracy = function(pred, truth) {
  mean(pred == truth)
}

#' f1
#' Calculate f1 based on comparing predicted and training labels
#' @param pred
#' @param truth
#'
#' @return f1
#' @export
#'
#' @examples

f1 = function(pred, truth) {
  p = precision(pred, truth)
  r = recall(pred, truth)
  2 * (p * r) / (p + r)
}

#' minimum_n_per_group
#' calculate the minimum n per group needed for a two-sample t-test
#' @param d
#' @param power
#'
#' @return the minimum n per group needed for a two-sample t-test
#' @export
#'
#' @examples

minimum_n_per_group = function(d, power = 0.8) {
  es = abs(d)
  qt = qt(1 - 0.5 * (1 - power), df = 2 * es^2)
  n = (2 * qt / d)^2
  return(ceiling(n))
}


#' r2
#' Calculate the r-squared statistic between predicted and ground truth continuous variables
#' @param pred
#' @param truth
#'
#' @return r-squared statistic between predicted and ground truth continuous variables
#' @export
#'
#' @examples

r2 = function(pred, truth) {
  ssr = sum((pred - truth)^2)
  sst = sum((truth - mean(truth))^2)
  return(1 - ssr / sst)
}


#' adj_R2
#' Calculate the adjusted r-squared statistic between predicted and ground truth continuous variables
#' @param pred
#' @param truth, n_p
#'
#' @return the adjusted r-squared statistic between predicted and ground truth continuous variables
#' @export
#'
#' @examples

adj_R2 = function(pred, truth, n_p) {
  n = length(truth)
  r2 = r2(pred, truth)
  adj_r2 = 1 - ((1 - r2) * (n - 1)) / (n - n_p - 1)
  return(adj_r2)
}





