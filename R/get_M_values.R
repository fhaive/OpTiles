#' Convert Beta Values to M-values
#'
#' This function transforms beta values (proportions in the range [0, 100]) to M-values using a logit2 (base-2 logit)
#' transformation. It applies a small \code{epsilon} correction to avoid issues with log(0) or division by zero.
#'
#' @param betavalue A \code{data.frame} or matrix containing beta values (typically in the range [0, 100]) for multiple samples.
#' @param epsilon A small numeric value added to avoid computing \code{log(0)} or \code{log(Inf)}. Default is \code{1e-2}.
#'
#' @return A \code{data.frame} with the same dimensions as \code{betavalue}, containing the M-value transformation of the input.
#'
#' @details The transformation is defined as:
#' \deqn{M = \log_2 \left(\frac{\beta}{1 - \beta} \right)}
#' where \eqn{\beta} values are constrained to the range [\code{epsilon}, 1 - \code{epsilon}] and scaled from [0, 100] to [0, 1].
#'
#' @examples
#' beta_matrix <- data.frame(
#'   sample1 = c(0, 50, 100),
#'   sample2 = c(10, 90, 50)
#' )
#' convert_beta_to_Mvalue(beta_matrix)
#'
#' @export


convert_beta_to_Mvalue <- function(betavalue, epsilon = 1e-2) {
  
  # Load required libraries
  library(data.table)
  library(dplyr)
  library(methylKit)
  library(ggplot2)
  library(tidyr)
  library(genomation)
  
  # Define logit transformation (base 2) for beta values
  logit2 <- function(beta) {
    beta <- beta/100 #because beta are in the range [0,100]
    #epsilon <- 1e-6  # Small value to prevent log(0) or division by zero
    beta <- pmax(pmin(beta, 1 - epsilon), epsilon)  # Keep beta values in valid range (0,1)
    log2(beta / (1 - beta))  # Apply logit transformation
  }
  
  # Convert beta values to M-values
  Mvalue <- sapply(betavalue, logit2)  # Apply logit function to each column (sample)
  # Add location column (preserve row names from original beta values)
  Mvalue <- as.data.frame(Mvalue)
  rownames(Mvalue) <- (rownames(betavalue))
  # Return transformed data frame
  return(Mvalue)
}
