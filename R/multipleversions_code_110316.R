## This file contains the functions for performing
## causal inference in randomized trials with control group
## heterogeneity
## Created: 11/3/16
## Modified: 11/3/16

#' Observed data quantities
#'
#' Calculates q_dyz = P(Y = y, D = d | Z = z)
#'
#' @param Y0 observed outcomes in Z = 0 arm
#' @param D0 observed treatment types in Z = 0 arm
#' @param Y1 observed outcomes in Z = 1 arm
#' @param D1 observed treatment types in Z = 1 arm
#' @export
get_observedprob <- function(Y0, D0, Y1, D1){

  # making sure that all inputs are numeric
  Y0 <- as.numeric(Y0)
  D0 <- as.numeric(D0)
  Y1 <- as.numeric(Y1)
  D1 <- as.numeric(D1)

  Y0 <- factor(Y0, levels = c(0,1))
  Y1 <- factor(Y1, levels = c(0,1))
  D0 <- factor(D0, levels = c(0,1,2))
  D1 <- factor(D1, levels = c(0,1))

  z0_prob_vals <- as.vector(table(D0, Y0))/length(Y0)

  names(z0_prob_vals) <- paste("q", outer(c(0,1,2),c(0,1), FUN = "paste", sep =""),0,sep = "" )

  z1_prob_vals <- as.vector(table(D1, Y1))/length(Y1)

  names(z1_prob_vals) <- paste("q",outer(c(0,1), c(0,1), FUN = "paste", sep = ""),1, sep = "" )

  good_probs <- c("q010","q100", "q110", "q200", "q210", "q001", "q011", "q111")

  z0_prob_vals <- z0_prob_vals[names(z0_prob_vals) %in% good_probs]
  z1_prob_vals <- z1_prob_vals[names(z1_prob_vals) %in% good_probs]

  # putting probabilities into specific order

  z0_prob_vals <- z0_prob_vals[c("q010","q100", "q110", "q200", "q210")]

  z1_prob_vals <- z1_prob_vals[c("q001", "q011", "q111")]

  c(z0_prob_vals, z1_prob_vals)
}


#' Asymptotic covariance matrix for estimatates of P(Y = y, D = d)
#'
#' @param Y vector of observed outcomes
#' @param D vector of observed treatment types
#' @param Z indicator for which treatment arm the data come from
#' @export
#' @return a matrix of variance-covariance estimates
#'

observedprob_vcov <- function(Y, D, Z){
  # This function computes the asymptotic covariance
  # matrix for estimates of P(Y = y, D = d)
  #
  # Args:
  #  Y: vector of observed outcomes
  #  D: vector of observed treatment types
  #
  # Returns: a covariance matrix
  #  with entries named "pdy" according to
  #  pdy = P(Y = y, D = d)

  prob_vals <- as.vector(table(D, Y))/length(Y)

  names(prob_vals) <- paste("q", outer(sort(levels(D)), sort(levels(Y)), FUN = "paste", sep =""),Z,sep = "" )


  # off-diagonals first
  cov_mat <- -prob_vals%o%prob_vals
  rownames(cov_mat) <- names(prob_vals)
  # putting variances on diagonal

  diag(cov_mat) <- prob_vals*(1-prob_vals)

  cov_mat
}


#' Gradient of upper and lower bounds on CACE
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D1 observed treatment types in the Z = 1 arm
#' @export
#' @return a 2 x 8 matrix of partial derivatives of the upper and lower bounds
#' with respect to the observed probabilites q_dyz
cace_bounds_grad <- function(Y0, D0, Y1, D1){

  args_list <- as.list(environment())


#   # constructing probabilities
#   z0_probs <- as.vector(table(D0, Y0)/length(Y0))
#   z1_probs <- as.vector(table(D1, Y1)/length(Y1))
#
#   # adding names
#   names(z0_probs) <- paste("q", outer(sort(unique(D0)), sort(unique(Y0)), FUN = "paste", sep = ""),
#                            "0", sep = "")
#   names(z1_probs) <- paste("q", outer(sort(unique(D1)), sort(unique(Y1)), FUN = 'paste', sep  = ""),
#                            "1", sep = "")
#
#   # removing probabiliites that aren't needed for the bounds
  good_probs <- c("q010","q100", "q110", "q200", "q210", "q001", "q011", "q111")
#
#   z0_probs <- z0_probs[names(z0_probs) %in% good_probs]
#   z1_probs <- z1_probs[names(z1_probs) %in% good_probs]

  # putting probabilities into specific order

  observed_probs <- do.call(get_observedprob, args_list)
  z0_probs <- observed_probs[1:5]

  z1_probs <- observed_probs[6:8]

  # constructing rows of derivative matrix
  w01 <- 1 - sum(z1_probs[c("q001", "q011")], z0_probs[c('q100', 'q110', 'q210', 'q200')])
  num <- z1_probs["q111"] - z0_probs["q110"] - z0_probs["q010"] + z1_probs["q011"]
  row2_num <- num - sum(z0_probs[c("q200","q210")])



  row1 <- c(-1/w01, num/w01^2, (-w01 + num)/w01^2, rep(num/w01^2, 3), (w01 + num)/w01^2, 1/w01)
  row2 <- c(-1/w01, row2_num/w01^2, rep((-w01 + row2_num)/w01^2, 3), row2_num/w01^2, (w01 + row2_num)/w01^2, 1/w01)


  results <- rbind(row1, row2)
  names(results) <- good_probs

  results

}


#' Gradient of bounds function
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D1 observed treatment types in the Z = 1 arm
#' @export
#' @return a 4 x 7 matrix of partial derivatives of the upper and lower bounds
#' with respect to the observed probabilites q_dyz
bounds_grad <- function(Y0, D0, Y1, D1){
  # Constructs the 4 x 7 matrix of partial derivatives
  # of the bounds with respect to the observed probabilities P(Y, D | Z)
  #
  # Args:
  #  Y0: observed outcomes in Z = 0 arm
  #  D0: observed treatment types in Z = 0 arm
  #  Y1: observed outcomes in Z = 1 arm
  #  D1: observed outcomes in Z = 1 arm
  #
  # Returns: a 4 x 7 matrix of partial derivatives

  # constructing probabilities
  D0 <- as.numeric(D0)
  D1 <- as.numeric(D1)
  Y0 <- as.numeric(Y0)
  Y1 <- as.numeric(Y1)

  z0_probs <- as.vector(table(D0, Y0)/length(Y0))
  z1_probs <- as.vector(table(D1, Y1)/length(Y1))

  # adding names
  names(z0_probs) <- paste("q", outer(sort(unique(D0)), sort(unique(Y0)), FUN = "paste", sep = ""),
                           "0", sep = "")
  names(z1_probs) <- paste("q", outer(sort(unique(D1)), sort(unique(Y1)), FUN = 'paste', sep  = ""),
                           "1", sep = "")

  # removing probabiliites that aren't needed for the bounds
  good_probs <- c("q010","q100", "q110", "q200", "q210", "q001", "q011", "q111")

  z0_probs <- z0_probs[names(z0_probs) %in% good_probs]
  z1_probs <- z1_probs[names(z1_probs) %in% good_probs]

  # putting probabilities into specific order

  z0_probs <- z0_probs[c("q100", "q110", "q200", "q210")]

  z1_probs <- z1_probs[c("q001", "q011", "q111")]

  # constructing rows of derivative matrix
  w01 <- 1 - sum(z1_probs[c("q001", "q011")], z0_probs[c('q100', 'q110', 'q210', 'q200')])
  q111_minus_q110 <- z1_probs["q111"] - z0_probs["q110"]
  w21 <- sum(z0_probs[c("q200","q210")])
  r1_common_entry <- q111_minus_q110/w01^2
  r2_common_entry <- w21/w01^2

  # gradient of w21/w01
  delta_w21overw01 <- c(rep(r2_common_entry, 2), rep(1/w01 + r2_common_entry, 2), rep(r2_common_entry, 2), 0)

  row1 <- c(r1_common_entry, r1_common_entry - 1/w01, rep(r1_common_entry, 4), 1/w01)
  row2 <- row1 - c(rep(r2_common_entry, 2), rep(1/w01 + r2_common_entry, 2), rep(r2_common_entry, 2), 0)
  row3 <- -delta_w21overw01*w01*q111_minus_q110/w21^2 + row1*w01/w21
  row4 <- row3 + delta_w21overw01*w01^2/w21^2

  results <- do.call(rbind, mget(paste("row", 1:4, sep = "")))
  names(results) <- good_probs

  results

}

#' Asymptotic covariance matrix of the upper and lower bounds estimators
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D1 observed treatment types in the Z = 1 arm
#' @export
#' @return a 4 x 4 estimated variance-covariance matrix
bounds_vcov <- function(Y0, D0, Y1, D1){
  # this function computes the asymptotic covariance matrix
  # of the upper and lower bounds for p_01 and p_21
  # based on the observed quantities P(Y = y, D = d | Z = z)
  #
  # Args:
  #  Y0: observed outcomes in the Z = 0 arm
  #  D0: observed treatment types in the Z = 1 arm
  #  Y1: observed outcomes in the Z = 1 arm
  #  D1: observed treatment types in the Z = 0 arm
  #
  # Returns: the 4 x 4 covariance matrix

  arg_list <- as.list(environment())

  # getting covariance matrix for Z = 0, 1 probabilities
  # then combining into one 10 x 10 matrix
  z1_covmat <- observedprob_vcov(Y = factor(Y1, levels = c(0,1)), D = factor(D1, levels = c(0,1)), Z = 1)

  z0_covmat <- observedprob_vcov(Y = factor(Y0, levels = c(0,1)), D = factor(D0, levels = c(0,1,2)), Z = 0)
  covmat <- rbind(z0_covmat, 0, 0, 0, 0)
  covmat <- cbind(covmat, 0, 0, 0, 0)
  covmat[7:10, 7:10] <- z1_covmat
  colnames(covmat) <- rownames(covmat) <- c(colnames(z0_covmat), colnames(z1_covmat))

  # trimming covariance matrix to only include the probabilities of interest
  good_probs <- c("q010","q100", "q110", "q200", "q210", "q001", "q011", "q111")
  covmat <- covmat[good_probs, good_probs]
  # getting gradient
  grad <- do.call(bounds_grad, arg_list)

  results <- grad%*%covmat%*%t(grad)

  colnames(results) <- rownames(results) <- c("p01_upper", "p01_lower", "p21_upper", "p21_lower")

  results

}


# Cn_finder <- function(n, delta, sigma, alpha, interval_vals = c(-10,10)){
#
#
#   Cn_func <- function(Cn){
#
#     pnorm(Cn + sqrt(n)*delta/max(sigma)) - pnorm(-Cn) - alpha
#
#   }
#
#   uniroot(Cn_func, interval = interval_vals)$root
# }

#' Bootstrapping upper and lower bounds
#'
#' Performs non-parameteric bootstrap of upper and lower bounds
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param D1 observed treatment types in the Z = 1 arm
#' @param n_boot number of bootstrap replications to perform
#' @export
#' @return an n_boot x 4 matrix of upper and lower bounds estimates
#' from the bootstrap samples
bounds_boot <- function(Y0, Y1, D0, D1, n_boot){
  # this function bootstraps the upper and lower bounds for p01 and p21
  #
  # Args:
  #  Y0: observed outcomes in the Z = 0 arm
  #  Y1: observed outcomes in the Z = 1 arm
  #  D0: observed treatment types in the Z = 0 arm
  #  D1: observed treatment types in the Z = 1 arm
  #  n_boot: the number of bootstrap replications to perform
  #
  # Returns:
  #  an n_boot x 4 matrix where each row is organized as
  #  (p01_lower, p01_upper, p21_lower, p21_upper)

  n0 <- length(Y0)
  n1 <- length(Y1)

  boot_Y0 <- rep(Y0, n_boot)
  boot_D0 <- rep(D0, n_boot)
  Z0_sample <- sample(rep(1:length(Y0), n_boot), replace = TRUE)

  Y0_mat <- matrix(boot_Y0[Z0_sample], nrow = n_boot, ncol = n0)
  D0_mat <- matrix(boot_D0[Z0_sample], nrow = n_boot, ncol = n0)

  boot_Y1 <- rep(Y1, n_boot)
  boot_D1 <- rep(D1, n_boot)
  Z1_sample <- sample(rep(1:n1, n_boot), replace = TRUE)

  Y1_mat <- matrix(boot_Y1[Z1_sample], nrow = n_boot, ncol = n1)
  D1_mat <- matrix(boot_D1[Z1_sample], nrow = n_boot, ncol = n1)

  Y1_D1_mat <- Y1_mat == 1 & D1_mat == 1
  Y0_D0_mat <- Y0_mat == 1 & D0_mat == 1

  q1_vec <- rowMeans(Y1_D1_mat)
  q0_vec <- rowMeans(Y0_D0_mat)

  omega_00_vec <- rowMeans(D1_mat == 0)
  omega_11_vec <- rowMeans(D0_mat == 1)
  omega_21_vec <- rowMeans(D0_mat == 2)

  omega_01_vec <- 1 - omega_00_vec - omega_11_vec - omega_21_vec

  p01_upper_vec <- (q1_vec - q0_vec)/omega_01_vec
  p01_lower_vec <- p01_upper_vec - omega_21_vec/omega_01_vec

  p21_upper_vec <- c(q1_vec - q0_vec)/omega_21_vec
  p21_lower_vec <- p21_upper_vec - omega_01_vec/omega_21_vec

  results <- cbind(p01_upper_vec, p01_lower_vec, p21_upper_vec, p21_lower_vec)

  colnames(results) <- c("p01_upper", "p01_lower", "p21_upper", "p21_lower")
  results

}

#' Find upper and lower bounds for P(Y(1) | D(0) = 0, D(1) = 1) and
#' P(Y(1) | D(0) = 2, D(1) = 1)
#'
#' @param Y0 observed outcomes in Z = 0 arm
#' @param Y1 observed outcomes in Z = 1 arm
#' @param D0 observed treatment types in D = 0 arm
#' @param D1 observed treatment types in D = 1 arm
#' @export
#' @return a 4-element vector of upper and lower bounds estimates
get_bounds <- function(Y0, Y1, D0, D1){
  # Finds the (possibly uninformative) bounds on
  # P(Y(1) = 1 | D(0) = 0, D(1) = 1) and P(Y(1) =1 | 2, 1)
  #
  # Args:
  #  Y0: observed outcomes in Z = 0 arm
  #  Y1: observed outcomes in Z = 1 arm
  #  D0: observed treatment type in Z = 0 arm
  #  D1: observed treatment type in Z = 1 arm
  #
  # Returns:
  #  a 4-element vector of upper and lower bounds estimates

  rhs <- mean(Y1 == 1 & D1 == 1) - mean(Y0 == 1 & D0 == 1)
  omega_00 <- mean(D1 == 0)
  omega_11 <- mean(D0 == 1)
  omega_21 <- mean(D0 == 2)
  omega_01 <- 1 - omega_11 - omega_21 - omega_00

  p01_lower <- (rhs - omega_21)/omega_01
  p01_upper <- rhs/omega_01

  p21_lower <- (rhs - omega_01)/omega_21
  p21_upper <- rhs/omega_21

  c(p01_upper = p01_upper, p01_lower =  p01_lower, p21_upper =  p21_upper,
    p21_lower = p21_lower)

}


#' Inference for Point Identified Quantities
#'
#' Finds point estimates, standard errors and confidence intervals
#' for point identified quantites in the multiple versions framework
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D1 observed outcomes in the Z = 1 arm
#'
#' @export
#' @return a matrix containing point estimates, standard errors
#' and confidence intervals
pointident_fit <- function(Y0, D0, Y1, D1, alpha_level){
  args_list <- list(Y0, D0, Y1, D1)
  observed_probs <- do.call(get_observedprob, args_list)
  observed_vcov_hat <- do.call(observedprob_vcov, args_list)

  p001 <- observed_probs["q011"]/(observed_probs["q011"] + observed_probs["q001"])
  p110 <- observed_probs["q110"]/(observed_probs["q110"] + observed_probs["q100"])
  p210 <- observed_probs["q210"]/(observed_probs["q210"] + observed_probs["q200"])

  omega_11 <- observed_probs["q100"] + observed_probs["q110"]
  omega_21 <- observed_probs["q200"] + observed_probs["q210"]
  omega_00 <- observed_probs["q011"] + observed_probs["q001"]
  omega_01 <- 1 - sum(omega_11, omega_21, omega_00)

  omega_vec <- c(omega_00, omega_01, omega_11, omega_21)

  omega_vcov <- -omega_vec%o%omega_vec
  diag(omega_vcov) <- omega_vec*(1 - omega_vec)


}

#' Causal inference with multiple versions of control
#'
#' Finds upper and lower bounds for the unidentified parameters
#' along with confidence intervals for the identified region
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param D1 observed treatment types in the Z = 1 arm
#' @param alpha_level the significance level at which confidence
#' intervals should be constructed
#' @export
#' @return a 4 x 4 matrix containing the bounds estimates, standard errors
#' and upper and lower confidence bounds
bounds_fit <- function(Y0, Y1, D0, D1, alpha_level){

  args_list <- list(Y0 = Y0, Y1 = Y1, D0 = D0, D1 = D1)

  bounds <- do.call(get_bounds, args_list)
  bounds_vcov_hat <- do.call(bounds_vcov, args_list)


  bounds_se <- sqrt(diag(bounds_vcov_hat)/N*2)

  bounds_ci <- bounds + bounds_se%o%qnorm(c(alpha_level/2, 1 - alpha_level/2))

  bound_names <- c("p01_upper", "p01_lower", "p21_upper", "p21_lower")
  bounds_fit <-  cbind(bounds, bounds_se, bounds_ci)

  colnames(bounds_fit) <- c("Estimate", "SE", "Lower Conf. Bound", "Upper Conf. Bound")

  rownames(bounds_fit) <- bound_names

  bounds_fit
}

#' Find upper and lower bounds for the CACE
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param D1 observed treatment types in the Z = 1 arm
#'
#' @export
#' @return a 2-element vector containing the upper and lower bounds for the CACE
get_cace <- function(Y0, Y1, D0, D1){

  args_list <- list(Y0 = Y0, Y1 = Y1, D0 = D0, D1 = D1)

  observed_probs <- do.call(get_observedprob, args_list)

  p01_bounds <- do.call(get_bounds, args_list)[1:2]

  p010 <- (observed_probs["q010"] - observed_probs["q011"])/(1 - sum(observed_probs[c("q200", "q210", "q100", "q110", "q001", "q011")]))

  output <- p01_bounds - p010
  names(output) <- c("cace_upper", "cace_lower")
  output

}

#' Find upper and lower bounds for the GCACE
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param D1 observed treatment types in the Z = 1 arm
#'
#' @export
#' @return a 2-element vector containing the upper and lower bounds for the GCACE
get_gcace <- function(Y0, Y1, D0, D1){

  args_list <- list(Y0 = Y0, Y1 = Y1, D0 = D0, D1 = D1)

  observed_probs <- do.call(get_observedprob, args_list)

  p21_bounds <- do.call(get_bounds, args_list)[3:4]

  p210 <- (observed_probs["q210"])/(observed_probs["q210"] + observed_probs["q200"])

  output <- p21_bounds - p210
  names(output) <- c("gcace_upper", "gcace_lower")
  output

}

#' Asymptotic covariance matrix for upper and lower bounds on CACE
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param D1 observed treatment types in the Z = 1 arm
#'
#' @export
#' @return a 2 x 2 covariance matrix for the upper and lower bounds on the CACE
cace_vcov <- function(Y0, Y1, D0, D1){

  arg_list <- as.list(environment())

  D0 <- as.numeric(D0)
  D1 <- as.numeric(D1)
  Y0 <- as.numeric(Y0)
  Y1 <- as.numeric(Y1)

  # getting covariance matrix for Z = 0, 1 probabilities
  # then combining into one 10 x 10 matrix

  z1_covmat <- observedprob_vcov(Y = factor(Y1, levels = c(0,1)), D = factor(D1, levels = c(0,1)), Z = 1)

  z0_covmat <- observedprob_vcov(Y = factor(Y0, levels = c(0,1)), D = factor(D0, levels = c(0,1,2)), Z = 0)
  covmat <- rbind(z0_covmat, 0, 0, 0, 0)
  covmat <- cbind(covmat, 0, 0, 0, 0)
  covmat[7:10, 7:10] <- z1_covmat
  colnames(covmat) <- rownames(covmat) <- c(colnames(z0_covmat), colnames(z1_covmat))

  # trimming covariance matrix to only include the probabilities of interest
  good_probs <- c("q010","q100", "q110", "q200", "q210", "q001", "q011", "q111")
  covmat <- covmat[good_probs, good_probs]
  # getting gradient
  grad <- do.call(cace_bounds_grad, arg_list)

  results <- grad%*%covmat%*%t(grad)

  colnames(results) <- rownames(results) <- c("cace_upper", "cace_lower")

  results



}

#' Testing model assumptions
#'
#' \code{model_test} checks whether the data violate the assumptions underlying the model,
#' similar to the assumptions tests for the instrumental variable model.
#'
#' @param Y0 observed outcomes in the Z = 0 arm
#' @param Y1 observed outcomes in the Z = 1 arm
#' @param D0 observed treatment types in the Z = 0 arm
#' @param D1 observed treatment types in the Z = 1 arm
#'
#' @export
#' @return a logical, TRUE indicates that the data are consistent with the assumptions,
#' FALSE indicates that the data violate at least one of the assumptions.
model_test <- function(Y0, Y1, D0, D1){

  test1 <- mean(Y0 == 1 & D0 == 0) - mean(Y1 == 1 & D1 == 0) >= 0
  test2 <- mean(Y1 == 1 & D1 == 1) - mean(Y0 == 1 & D0 == 1) >= 0

  test1 & test2

}
